classdef gas_cylinder_model < addition_gas_cylinder_parallel_valve
%gas_cylinder_model  三参数构造别名，兼容脚本调用：
%  cyl = gas_cylinder_model(gas_initial_conditions, orifice, [valve1, valve2]);
%内部委托 addition_gas_cylinder_parallel_valve（孔板 + 两只并联电磁阀 + 壁面换热）。
%initial_conditions 压力/温度字段可用 PO_init/TO_init 或 P0_init/T0_init。

methods
    function obj = gas_cylinder_model(initial_conditions, orifice, valves)
        if nargin ~= 3
            error('gas_cylinder_model:BadNargin', ...
                ['用法: gas_cylinder_model(initial_conditions, orifice, [valve1, valve2])。', ...
                '需要恰好 3 个输入参数。']);
        end
        obj@addition_gas_cylinder_parallel_valve(initial_conditions, orifice, valves);
    end
end
end
