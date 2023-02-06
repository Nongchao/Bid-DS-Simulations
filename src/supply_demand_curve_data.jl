# run through line 91 of "bid_model_driver.jl" to get the system
using Plots

rep_time = 13
base_power = get_base_power(sys)
total_load = 0.0
renewable_dispatch = 0.0
renewable_fix = 0.0
hydro_dispatch = 0.0
hydro_reservoir = 0.0

for l in get_components(PowerLoad, sys)
    ts = values(get_time_series_array(SingleTimeSeries, l, "max_active_power"))
    total_load = total_load + ts[rep_time] * base_power
end

for g in get_components(RenewableDispatch, sys)
    ts = values(get_time_series_array(SingleTimeSeries, g, "max_active_power"))
    renewable_dispatch = renewable_dispatch + ts[rep_time] * base_power
end

for g in get_components(RenewableFix, sys)
    ts = values(get_time_series_array(SingleTimeSeries, g, "max_active_power"))
    renewable_fix = renewable_fix + ts[rep_time] * base_power
end

for g in get_components(HydroDispatch, sys)
    ts = values(get_time_series_array(SingleTimeSeries, g, "max_active_power"))
    hydro_dispatch = hydro_dispatch + ts[rep_time] * base_power
end

for g in get_components(HydroEnergyReservoir, sys)
    ts = values(get_time_series_array(SingleTimeSeries, g, "max_active_power"))
    hydro_reservoir = hydro_reservoir + ts[rep_time] * base_power
end

variable_cost_go = []
for g in get_components(ThermalStandard, sys)
    variable_cost_psy = get_variable(get_operation_cost(g)).cost

    for b in 1:length(variable_cost_psy)
        if b == 1
            if length(variable_cost_psy[b]) == 1
                # for this sys, this is the case of the "<1|2|3>14_SYNC_COND_1" devices with zero cost
                block = [get_variable(get_operation_cost(g)).cost, get_max_active_power(g) * base_power]
            else
                block = [variable_cost_psy[b][1]/variable_cost_psy[b][2], variable_cost_psy[b][2]]
            end
        else
            block = [(variable_cost_psy[b][1]-variable_cost_psy[b-1][1])/(variable_cost_psy[b][2]-variable_cost_psy[b-1][2]), (variable_cost_psy[b][2]-variable_cost_psy[b-1][2])]
        end
        push!(variable_cost_go, block)
    end
end

push!(variable_cost_go, [0.0, renewable_dispatch + renewable_fix + hydro_dispatch + hydro_reservoir])

supply_stack_df = DataFrame([[], []], ["var_cost", "quantity"])

for i in 1:length(variable_cost_go)
    push!(supply_stack_df, [variable_cost_go[i][1], variable_cost_go[i][2]])
end

sort!(supply_stack_df, [:var_cost])