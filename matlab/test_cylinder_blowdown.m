% test_cylinder_blowdown  物理正确性检查：壁面换热气瓶放气至 ~0.5 MPa。
%
% 用法（在仓库根目录或本目录）:
%   addpath(fullfile(pwd, 'matlab'));   % 或 addpath('<repo>/matlab')
%   test_cylinder_blowdown
%
% 通过准则（有 REFPROP 时）:
%   - 放气至约 0.5 MPa 后气瓶温度 T_final > T_min_pass（默认 200 K）
%   - T_final 明显高于同压比等熵估计（理想气体或 REFPROP S=const）
% 无 REFPROP 时:
%   - 仍做构造 + networkFlow 冒烟；跳过 mass_flow 状态更新并打印说明

function test_cylinder_blowdown()
    thisDir = fileparts(mfilename('fullpath'));
    if isempty(thisDir)
        thisDir = pwd;
    end
    addpath(thisDir);

    P0 = 35e6;          % Pa
    T0 = 300;           % K
    V = 0.18;           % m^3
    P_stop = 0.5e6;     % Pa — 与模型内截止一致
    tankP = 0.1e6;      % Pa 背压
    dt = 0.01;          % s 初始步长
    max_steps = 200000;
    T_min_pass = 200;   % K — 不得落到“几十 K”
    log_every = 50;

    gas_initial_conditions = struct();
    gas_initial_conditions.P0_init = P0;
    gas_initial_conditions.T0_init = T0;
    gas_initial_conditions.V = V;

    orifice = struct('Cd', 1, 'A', pi * (0.005)^2);   % GLQ1035（按用户：半径 5 mm）
    valve1 = struct('Cd', 0.9, 'A', pi * (0.001)^2); % DC22F2035
    valve2 = struct('Cd', 0.9, 'A', pi * (0.001)^2);

    fprintf('=== test_cylinder_blowdown ===\n');
    fprintf('Init: P0=%.3g Pa, T0=%.1f K, V=%.3f m^3\n', P0, T0, V);

    cyl = gas_cylinder_model(gas_initial_conditions, orifice, [valve1, valve2]);
    fprintf('Construct OK: class=%s, V=%.3f, h_conv=%.1f W/(m^2·K)\n', ...
        class(cyl), cyl.V, cyl.h_conv);

    % networkFlow 冒烟（理想气体孔板/阀，不需 REFPROP）
    flow0 = cyl.networkFlow(P0, T0, tankP, true, true);
    fprintf('networkFlow smoke: m_dot(P0,T0)=%.4g kg/s (valves open)\n', flow0);
    assert(flow0 > 0, 'networkFlow should be positive with both valves open');
    flow_closed = cyl.networkFlow(P0, T0, tankP, false, false);
    assert(flow_closed == 0, 'networkFlow should be 0 with both valves closed');

    hasRefprop = ~isempty(which('refpropm'));
    if ~hasRefprop
        fprintf('\nSKIP mass_flow blowdown: refpropm (REFPROP) not on MATLAB path.\n');
        fprintf('  Construct + networkFlow smoke PASSED.\n');
        fprintf('  To run the heat-transfer state-update check, install REFPROP and add refpropm.\n');
        return
    end

    % --- 等熵对照（同压比，不推进质量）---
    gamma = cyl.gamma;
    T_isen_ideal = T0 * (P_stop / P0) ^ ((gamma - 1) / gamma);
    T_isen_ref = NaN;
    try
        S0 = addition_gas_cylinder_parallel_valve.heliumEntropy(T0, P0);
        T_isen_ref = fzero(@(T) refpropm('S', 'T', T, 'P', P_stop / 1000, 'Helium') - S0, ...
            max(T_isen_ideal, 40));
    catch ME
        fprintf('Note: REFPROP isentropic estimate failed (%s); using ideal-gas only.\n', ME.message);
    end
    fprintf('Isentropic estimate to %.2f MPa: ideal-gas T≈%.1f K', P_stop / 1e6, T_isen_ideal);
    if ~isnan(T_isen_ref)
        fprintf(', REFPROP S=const T≈%.1f K', T_isen_ref);
    end
    fprintf('\n');

    % --- 换热模型放气积分 ---
    conditions = gas_initial_conditions;
    t = 0;
    n = 0;
    log_t = t;
    log_P = P0;
    log_T = T0;
    log_m = flow0;

    fprintf('Stepping mass_flow until P <= %.2f MPa (dt=%.3g s)...\n', P_stop / 1e6, dt);
    while n < max_steps
        P = localPressure(conditions);
        if P <= P_stop
            break
        end

        % 限制单步质量取出，避免抽空
        rho = addition_gas_cylinder_parallel_valve.heliumDensity(localTemp(conditions), P);
        m_now = rho * V;
        flow_est = cyl.networkFlow(P, localTemp(conditions), tankP, true, true);
        if flow_est <= 0
            fprintf('Stop: networkFlow returned <= 0 at P=%.3g Pa\n', P);
            break
        end
        dt_use = dt;
        if flow_est * dt_use > 0.02 * m_now
            dt_use = 0.02 * m_now / flow_est;
        end

        [result, conditions] = cyl.mass_flow(conditions, tankP, dt_use, true, true);
        t = t + dt_use;
        n = n + 1;
        flow = result(1);
        T_exit = result(2);           % 节流出口温度（非瓶内主体温度）
        T_cyl = localTemp(conditions); % 瓶内温度在 new_conditions
        P_now = localPressure(conditions);

        if mod(n, log_every) == 0 || P_now <= P_stop
            fprintf('  t=%8.3f s  P=%8.3f MPa  T_cyl=%7.2f K  T_exit=%7.2f K  m_dot=%8.4g kg/s\n', ...
                t, P_now / 1e6, T_cyl, T_exit, flow);
        end

        log_t(end + 1, 1) = t; %#ok<AGROW>
        log_P(end + 1, 1) = P_now; %#ok<AGROW>
        log_T(end + 1, 1) = T_cyl; %#ok<AGROW>
        log_m(end + 1, 1) = flow; %#ok<AGROW>

        if P_now <= P_stop
            break
        end
        if flow <= 0 && P_now > P_stop
            warning('test_cylinder_blowdown:ZeroFlow', ...
                'Flow became zero before P_stop; last P=%.3g Pa', P_now);
            break
        end
    end

    if n >= max_steps
        error('test_cylinder_blowdown:MaxSteps', ...
            'Exceeded max_steps=%d before reaching %.2f MPa', max_steps, P_stop / 1e6);
    end

    P_final = localPressure(conditions);
    T_final = localTemp(conditions);

    fprintf('\n--- Results ---\n');
    fprintf('Steps=%d, t_final=%.3f s\n', n, t);
    fprintf('Final state: P=%.4f MPa, T=%.2f K\n', P_final / 1e6, T_final);
    fprintf('Pass floor: T > %.0f K (not tens of kelvin)\n', T_min_pass);
    fprintf('Compare: heat-transfer T=%.1f K vs isentropic ideal≈%.1f K', T_final, T_isen_ideal);
    if ~isnan(T_isen_ref)
        fprintf(' / REFPROP≈%.1f K', T_isen_ref);
    end
    fprintf('\n');

    assignin('base', 'blowdown_log', struct( ...
        't', log_t, 'P', log_P, 'T', log_T, 'flow', log_m, ...
        'P_final', P_final, 'T_final', T_final, ...
        'T_isen_ideal', T_isen_ideal, 'T_isen_refprop', T_isen_ref));

    assert(T_final > T_min_pass, ...
        'FAIL: T_final=%.2f K <= %.0f K (looks like near-isentropic / unphysical cooling)', ...
        T_final, T_min_pass);

    T_isen_ref_for_cmp = T_isen_ideal;
    if ~isnan(T_isen_ref)
        T_isen_ref_for_cmp = T_isen_ref;
    end
    if T_final <= T_isen_ref_for_cmp + 30
        warning('test_cylinder_blowdown:WeakHeatTransfer', ...
            ['T_final (%.1f K) is only slightly above isentropic (%.1f K). ', ...
            'Defaults should keep gas much warmer; check h_conv / A_wall.'], ...
            T_final, T_isen_ref_for_cmp);
    end

    fprintf('PASS: heat-transfer blowdown T_final=%.1f K > %.0f K (isentropic would be ~%.0f K).\n', ...
        T_final, T_min_pass, T_isen_ref_for_cmp);
end

function P = localPressure(conditions)
    if isfield(conditions, 'PO_init') && ~isempty(conditions.PO_init)
        P = conditions.PO_init;
    else
        P = conditions.P0_init;
    end
end

function T = localTemp(conditions)
    if isfield(conditions, 'TO_init') && ~isempty(conditions.TO_init)
        T = conditions.TO_init;
    else
        T = conditions.T0_init;
    end
end
