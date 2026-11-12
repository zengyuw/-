classdef addition_gas_cylinder_parallel_valve < handle
%后续考虑将gamma与R用refpropm调用
%气瓶放气状态更新采用开系能量衡算 + 壁面/环境换热，避免纯等熵导致的过冷。
properties
    gamma = 1.67
    R = 2077
    orifice
    valves
    V
    P_junction = []
    branch_flows = [0, 0]
    valve_open = [false, false]
    P_critical_ratio
    % --- 壁面换热参数（可在 initial_conditions 中覆盖）---
    % h_conv [W/(m^2·K)]  气–内壁对流换热系数
    h_conv = 100
    % A_wall [m^2]  内壁换热面积；空则按气瓶容积估算
    A_wall = []
    % T_amb [K]  环境温度
    T_amb = 293.15
    % T_wall [K]  当前壁温（无壁面热容时保持为 T_amb）
    T_wall = 293.15
    % m_wall [kg]  壁面热质量；0 表示等温壁（无限热容/与环境充分耦合）
    m_wall = 0
    % c_wall [J/(kg·K)]  壁面比热（钢约 500）
    c_wall = 500
    % h_amb [W/(m^2·K)]  外壁–环境对流换热系数
    h_amb = 10
    % A_ext [m^2]  外壁换热面积；空则取与 A_wall 相同
    A_ext = []
end

methods
    function obj = addition_gas_cylinder_parallel_valve(initial_conditions, orifice, valves)
        if iscell(valves)
            valves = [valves{:}];
        end
        if ~isstruct(valves) || numel(valves) ~= 2
            error('addition_gas_cylinder_parallel_valve:BadValves', ...
                '需要两只并联电磁阀，valves = [电磁阀1，电磁阀2]');
        end
        if ~isfield(orifice, 'Cd') || ~isfield(orifice, 'A') ...
                || ~isfield(valves, 'Cd') || ~isfield(valves, 'A')
            error('addition_gas_cylinder_parallel_valve:BadElement', ...
                '孔板和电磁阀都需要字段 Cd 和 A');
        end
        obj.orifice = orifice;
        obj.valves = valves;
        obj.P_critical_ratio = (2 / (obj.gamma + 1)) ^ (obj.gamma / (obj.gamma - 1));
        obj.V = initial_conditions.V;
        obj.applyHeatTransferDefaults(initial_conditions);
    end

    function [result, new_conditions] = mass_flow(obj, conditions, tankP, time_step, valve1_open, valve2_open)
        %按电磁阀开闭信号计算总流量，并更新气瓶状态。
        %两阀都关闭，或气瓶压力低于0.5 MPa时不放气，状态保持不变。
        %状态更新：开系能量衡算 d(mu)/dt = -h_exit*m_dot + Q_dot，
        %其中 Q_dot = h_conv*A_wall*(T_wall - T_gas)，可选壁面热容。
        PO = conditions.PO_init;
        TO = conditions.TO_init;
        valve1_open = addition_gas_cylinder_parallel_valve.asOpen(valve1_open);
        valve2_open = addition_gas_cylinder_parallel_valve.asOpen(valve2_open);
        if PO < 0.5e6 || (~valve1_open && ~valve2_open)
            obj.P_junction = PO;
            obj.branch_flows = [0, 0];
            obj.valve_open = [valve1_open, valve2_open];
            result = [0, TO];
            new_conditions = conditions;
            return
        end
        flow = obj.networkFlow(PO, TO, tankP, valve1_open, valve2_open);
        delta_m = flow * time_step; %补压气体质量减少
        rho = addition_gas_cylinder_parallel_valve.heliumDensity(TO, PO);
        m0 = rho * obj.V;
        m_new = m0 - delta_m;
        if m_new <= 0
            error('addition_gas_cylinder_parallel_valve:EmptyCylinder', ...
                '本步放气质量超过气瓶存气，请减小 time_step 或检查流量。');
        end
        rho_new = m_new / obj.V;

        T_wall = obj.currentWallTemperature(conditions);
        A_wall = obj.effectiveWallArea();
        A_ext = obj.effectiveExtArea(A_wall);

        % 开系能量衡算（显式欧拉，出口焓取瓶内状态焓）
        u0 = addition_gas_cylinder_parallel_valve.heliumInternalEnergy(TO, PO);
        h_exit = addition_gas_cylinder_parallel_valve.heliumEnthalpy(TO, PO);
        Q_dot = obj.h_conv * A_wall * (T_wall - TO);
        Q = Q_dot * time_step;
        U_new = m0 * u0 - delta_m * h_exit + Q;
        u_new = U_new / m_new;

        TO_guess = max(TO, 50);
        TO_new = fzero(@(T) refpropm('U', 'T', T, 'D', rho_new, 'Helium') - u_new, TO_guess);
        PO_new = refpropm('P', 'T', TO_new, 'D', rho_new, 'Helium') * 1000;

        T_wall_new = obj.updateWallTemperature(T_wall, TO, time_step, A_wall, A_ext);
        obj.T_wall = T_wall_new;

        new_conditions = struct( ...
            'PO_init', max(PO_new, 0), ...
            'TO_init', max(TO_new, 0), ...
            'V', obj.V, ...
            'T_wall', T_wall_new);
        result = [max(flow, 0), TO_new];
    end

    function flow = networkFlow(obj, PO, TO, Pb, valve1_open, valve2_open)
        %求节点压力，使孔板流量等于两只电磁阀流量之和。返回离开气瓶的流量。
        valve1_open = addition_gas_cylinder_parallel_valve.asOpen(valve1_open);
        valve2_open = addition_gas_cylinder_parallel_valve.asOpen(valve2_open);
        obj.P_junction = PO;
        obj.branch_flows = [0, 0];
        obj.valve_open = [valve1_open, valve2_open];
        if (~valve1_open && ~valve2_open) || PO <= Pb || TO <= 0 || (PO - Pb) < 1
            flow = 0;
            return
        end
        span = PO - Pb;
        lo = Pb + 1e-8 * span;
        hi = PO - 1e-8 * span;
        residual = @(p_mid) obj.junctionResidual(p_mid, PO, TO, Pb, valve1_open, valve2_open);
        fLo = residual(lo);
        fHi = residual(hi);
        if fLo == 0
            p_mid = lo;
        elseif fHi == 0
            p_mid = hi;
        elseif fLo * fHi > 0
            if abs(fLo) <= abs(fHi)
                p_mid = lo;
            else
                p_mid = hi;
            end
        else
            p_mid = fzero(residual, [lo, hi]);
        end
        [m_front, t_mid] = obj.elementFlow(PO, TO, p_mid, obj.orifice);
        [m1, ~] = obj.valveFlow(p_mid, t_mid, Pb, obj.valves(1), valve1_open);
        [m2, ~] = obj.valveFlow(p_mid, t_mid, Pb, obj.valves(2), valve2_open);
        obj.P_junction = p_mid;
        obj.branch_flows = [m1, m2];
        flow = max(m_front, 0);
    end
end

methods (Access = private)
    function applyHeatTransferDefaults(obj, ic)
        if isfield(ic, 'h_conv') && ~isempty(ic.h_conv)
            obj.h_conv = ic.h_conv;
        end
        if isfield(ic, 'A_wall') && ~isempty(ic.A_wall)
            obj.A_wall = ic.A_wall;
        end
        if isfield(ic, 'T_amb') && ~isempty(ic.T_amb)
            obj.T_amb = ic.T_amb;
        end
        if isfield(ic, 'T_wall') && ~isempty(ic.T_wall)
            obj.T_wall = ic.T_wall;
        else
            obj.T_wall = obj.T_amb;
        end
        if isfield(ic, 'm_wall') && ~isempty(ic.m_wall)
            obj.m_wall = ic.m_wall;
        end
        if isfield(ic, 'c_wall') && ~isempty(ic.c_wall)
            obj.c_wall = ic.c_wall;
        end
        if isfield(ic, 'h_amb') && ~isempty(ic.h_amb)
            obj.h_amb = ic.h_amb;
        end
        if isfield(ic, 'A_ext') && ~isempty(ic.A_ext)
            obj.A_ext = ic.A_ext;
        end
    end

    function T_wall = currentWallTemperature(obj, conditions)
        if isfield(conditions, 'T_wall') && ~isempty(conditions.T_wall)
            T_wall = conditions.T_wall;
        else
            T_wall = obj.T_wall;
        end
    end

    function A_wall = effectiveWallArea(obj)
        if isempty(obj.A_wall) || obj.A_wall <= 0
            A_wall = addition_gas_cylinder_parallel_valve.estimateWallArea(obj.V);
        else
            A_wall = obj.A_wall;
        end
    end

    function A_ext = effectiveExtArea(obj, A_wall)
        if isempty(obj.A_ext) || obj.A_ext <= 0
            A_ext = A_wall;
        else
            A_ext = obj.A_ext;
        end
    end

    function T_wall_new = updateWallTemperature(obj, T_wall, T_gas, time_step, A_wall, A_ext)
        % m_wall==0：等温壁，壁温维持环境温度（强热耦合极限）
        if obj.m_wall <= 0 || obj.c_wall <= 0
            T_wall_new = obj.T_amb;
            return
        end
        Q_gas = obj.h_conv * A_wall * (T_wall - T_gas) * time_step;
        Q_amb = obj.h_amb * A_ext * (obj.T_amb - T_wall) * time_step;
        T_wall_new = T_wall + (-Q_gas + Q_amb) / (obj.m_wall * obj.c_wall);
    end

    function [flow, T_out] = elementFlow(obj, P_in, T_in, P_out, element)
        if P_in <= 0 || T_in <= 0 || P_out >= P_in
            flow = 0;
            T_out = T_in;
            return
        end
        Cd = element.Cd;
        A = element.A;
        pr = P_out / P_in;
        if pr <= obj.P_critical_ratio
            term = (2 / (obj.gamma + 1)) ^ ((obj.gamma + 1) / (obj.gamma - 1));
            flow = Cd * A * P_in * sqrt(obj.gamma / (obj.R * T_in)) * sqrt(term);
            T_out = T_in * (2 / (obj.gamma + 1));
        else
            term = pr ^ (2 / obj.gamma) - pr ^ ((obj.gamma + 1) / obj.gamma);
            flow = Cd * A * P_in * sqrt(2 * obj.gamma / ((obj.gamma - 1) * obj.R * T_in) * term);
            T_out = T_in * pr ^ ((obj.gamma - 1) / obj.gamma);
        end
        flow = max(flow, 0);
    end

    function [flow, T_out] = valveFlow(obj, P_in, T_in, P_out, valve, is_open)
        if ~is_open
            flow = 0;
            T_out = T_in;
            return
        end
        [flow, T_out] = obj.elementFlow(P_in, T_in, P_out, valve);
    end

    function residual = junctionResidual(obj, p_mid, P0, T0, Pb, valve1_open, valve2_open)
        [m_front, t_mid] = obj.elementFlow(P0, T0, p_mid, obj.orifice);
        [m1, ~] = obj.valveFlow(p_mid, t_mid, Pb, obj.valves(1), valve1_open);
        [m2, ~] = obj.valveFlow(p_mid, t_mid, Pb, obj.valves(2), valve2_open);
        residual = m_front - (m1 + m2);
    end
end

methods (Static)
    function rho = heliumDensity(T, P)
        % refpropm 压力单位为 kPa，这里的 P 为 Pa。
        addition_gas_cylinder_parallel_valve.requireRefprop();
        rho = double(refpropm('D', 'T', T, 'P', P / 1000, 'Helium'));
    end

    function S = heliumEntropy(T, P)
        addition_gas_cylinder_parallel_valve.requireRefprop();
        S = double(refpropm('S', 'T', T, 'P', P / 1000, 'Helium'));
    end

    function u = heliumInternalEnergy(T, P)
        % 比内能 [J/kg]；P 为 Pa。
        addition_gas_cylinder_parallel_valve.requireRefprop();
        u = double(refpropm('U', 'T', T, 'P', P / 1000, 'Helium'));
    end

    function h = heliumEnthalpy(T, P)
        % 比焓 [J/kg]；P 为 Pa。
        addition_gas_cylinder_parallel_valve.requireRefprop();
        h = double(refpropm('H', 'T', T, 'P', P / 1000, 'Helium'));
    end

    function A = estimateWallArea(V)
        % 按长径比 L/D ≈ 4 的圆筒估算内壁面积 [m^2]。
        % V = π r^2 L，L = 8 r → r = (V/(8π))^(1/3)，A = 18 π r^2。
        if ~(isscalar(V) && isnumeric(V) && isreal(V) && V > 0)
            error('addition_gas_cylinder_parallel_valve:BadVolume', ...
                '气瓶容积 V 须为正实数 [m^3]');
        end
        r = (V / (8 * pi)) ^ (1 / 3);
        A = 18 * pi * r ^ 2;
    end
end

methods (Static, Access = private)
    function requireRefprop()
        if isempty(which('refpropm'))
            error('addition_gas_cylinder_parallel_valve:NoRefprop', ...
                ['气瓶状态更新需要 refpropm (REFPROP 的 MATLAB 接口)。', ...
                '孔板和电磁阀流量可调用 networkFlow，不需要物性库。']);
        end
    end

    function flag = asOpen(signal)
        if ~(isscalar(signal) && (islogical(signal) || isnumeric(signal))) || ~isreal(signal)
            error('addition_gas_cylinder_parallel_valve:BadSignal', ...
                '电磁阀信号应为 0/1 或 true/false');
        end
        flag = signal ~= 0;
    end
end
end
