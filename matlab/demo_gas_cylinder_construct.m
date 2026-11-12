% demo_gas_cylinder_construct  与用户脚本一致的气瓶构造示例。
% 将本目录加入路径后运行本脚本即可构造 cyl（不调用 REFPROP）。
%
% 注意：orifice.A = pi*(0.005)^2 按用户脚本（注释为半径 5 mm）。

gas_initial_conditions = struct();
gas_initial_conditions.P0_init = 35e6;  % 气瓶初始压力 Pa
gas_initial_conditions.T0_init = 300;   % 气瓶初始温度 K
gas_initial_conditions.V = 0.18;        % 气瓶容积 m^3
conditions = gas_initial_conditions;

orifice = struct('Cd', 1, 'A', pi * (0.005)^2);   % GLQ1035
valve1 = struct('Cd', 0.9, 'A', pi * (0.001)^2); % DC22F2035
valve2 = struct('Cd', 0.9, 'A', pi * (0.001)^2); % DC22F2035

cyl = gas_cylinder_model(gas_initial_conditions, orifice, [valve1, valve2]);

fprintf('cyl class: %s\n', class(cyl));
fprintf('V = %.3f m^3, orifice.A = %.3e m^2\n', cyl.V, cyl.orifice.A);
disp('Construction OK. Example mass_flow call (needs REFPROP):');
disp('[result, conditions] = cyl.mass_flow(conditions, tankP, dt, true, true);');
