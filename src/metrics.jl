function results_output(
    res::PSI.SimulationProblemResults{PSI.DecisionModelSimulationResults}, 
    sim::Simulation,
    PTDF_matrix::PTDF{Tuple{Vector{String}, Vector{Int64}}, Tuple{Dict{String, Int64}, Dict{Int64, Int64}}}, 
    emission_type::Vector{String}, 
    source_gen_data::DataFrame,
    result_file_path::String,
    start_date::Dates.DateTime,
    load_output_model::String)
    
    duals_line = PSI.read_realized_dual(res, "NetworkFlowConstraint__Line")
    duals_taptransformer = PSI.read_realized_dual(res, "NetworkFlowConstraint__TapTransformer")
    λ = PSI.read_realized_dual(res, "CopperPlateBalanceConstraint__System")[:, 2]
    flow_duals = DataFrames.outerjoin(
        values(duals_line),
        values(duals_taptransformer),
        on = :DateTime,
    )
    μ = Matrix(flow_duals[:, PTDF_matrix.axes[1]])

    sys = PSI.get_system(res)
    base_power = get_base_power(get_system(res))
    LMP = flow_duals[:, [:DateTime]]
    for bus in get_components(Bus, sys)
        LMP[:, get_name(bus)] = (λ .+ μ * PTDF_matrix[:, get_number(bus)]) ./ base_power
    end

    CSV.write(joinpath(result_file_path, "lmp.csv"), LMP)

    # reserve prices
    # p_Reg_Up = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveUp__Reg_Up")[:, 2] ./ base_power
    # p_Reg_Down = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveDown__Reg_Down")[:, 2] ./ base_power
    # p_Flex_Up = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveUp__Flex_Up")[:, 2] ./ base_power
    # p_Flex_Down = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveDown__Flex_Down")[:, 2] ./ base_power
    # p_Spin_Up_R1 = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveUp__Spin_Up_R1")[:, 2] ./ base_power
    # p_Spin_Up_R2 = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveUp__Spin_Up_R2")[:, 2] ./ base_power
    # p_Spin_Up_R3 = PSI.read_realized_dual(res, "RequirementConstraint__VariableReserve__ReserveUp__Spin_Up_R3")[:, 2] ./ base_power

    variable_values = PSI.read_realized_variables(res)
    optimizer_stats = PSI.read_optimizer_stats(res)
    CSV.write(joinpath(result_file_path, "optimizer_stats.csv"), optimizer_stats)

    thermal_gen = variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "thermal_generation.csv"), thermal_gen)
    renewable_dispatch = variable_values["ActivePowerVariable__RenewableDispatch"][!, names(variable_values["ActivePowerVariable__RenewableDispatch"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "renewable_dispatch.csv"), renewable_dispatch)
    hydro_dispatch = variable_values["ActivePowerVariable__HydroDispatch"][!, names(variable_values["ActivePowerVariable__HydroDispatch"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "hydro_dispatch.csv"), hydro_dispatch)
    hydro_reservoir = variable_values["ActivePowerVariable__HydroEnergyReservoir"][!, names(variable_values["ActivePowerVariable__HydroEnergyReservoir"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "hydro_reservoir.csv"), hydro_reservoir)
    if load_output_model == "bid_in_load"
        flex_load = variable_values["ActivePowerVariable__FlexibleLoad"][!, names(variable_values["ActivePowerVariable__FlexibleLoad"]) .!= "DateTime"]
        CSV.write(joinpath(result_file_path, "flex_load_dispatch.csv"), flex_load)
    elseif load_output_model == "shiftable_load"
        shift_load = variable_values["ActivePowerVariable__ShiftableLoad"][!, names(variable_values["ActivePowerVariable__ShiftableLoad"]) .!= "DateTime"]
        CSV.write(joinpath(result_file_path, "shift_load.csv"), shift_load)
    end
    line_flow = variable_values["FlowActivePowerVariable__Line"][!, names(variable_values["FlowActivePowerVariable__Line"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "line_flow.csv"), line_flow)
    transformer_flow = variable_values["FlowActivePowerVariable__TapTransformer"][!, names(variable_values["FlowActivePowerVariable__TapTransformer"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "transformer_flow.csv"), transformer_flow)
    thermal_gen_on = variable_values["OnVariable__ThermalStandard"][!, names(variable_values["OnVariable__ThermalStandard"]) .!= "DateTime"]
    CSV.write(joinpath(result_file_path, "thermal_gen_on.csv"), thermal_gen_on)

    time_step = size(thermal_gen)[1]

    ################### total reserve payment and cost ########################
    # q_Reg_Up = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Reg_Up"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Reg_Up"]) .!= "DateTime"]
    # q_Reg_Down = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveDown__Reg_Down"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveDown__Reg_Down"]) .!= "DateTime"]
    # q_Flex_Up = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Flex_Up"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Flex_Up"]) .!= "DateTime"]
    # q_Flex_Down = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveDown__Flex_Down"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveDown__Flex_Down"]) .!= "DateTime"]
    # q_Spin_Up_R1 = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Spin_Up_R1"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Spin_Up_R1"]) .!= "DateTime"]
    # q_Spin_Up_R2 = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Spin_Up_R2"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Spin_Up_R2"]) .!= "DateTime"]
    # q_Spin_Up_R3 = variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Spin_Up_R3"][!, names(variable_values["ActivePowerReserveVariable__VariableReserve__ReserveUp__Spin_Up_R3"]) .!= "DateTime"]

    # reserve_payment = sum(sum(eachcol(p_Reg_Up .* q_Reg_Up))) + 
    #     sum(sum(eachcol(p_Reg_Down .* q_Reg_Down))) + 
    #     sum(sum(eachcol(p_Flex_Up .* q_Flex_Up))) + 
    #     sum(sum(eachcol(p_Flex_Down .* q_Flex_Down))) + 
    #     sum(sum(eachcol(p_Spin_Up_R1 .* q_Spin_Up_R1))) + 
    #     sum(sum(eachcol(p_Spin_Up_R2 .* q_Spin_Up_R2))) + 
    #     sum(sum(eachcol(p_Spin_Up_R3 .* q_Spin_Up_R3)))

    # reserve_cost = sum(sum(eachcol(0.0001 .* q_Reg_Up))) + 
    #     sum(sum(eachcol(0.0001 .* q_Reg_Down))) + 
    #     sum(sum(eachcol(0.0001 .* q_Flex_Up))) + 
    #     sum(sum(eachcol(0.0001 .* q_Flex_Down))) + 
    #     sum(sum(eachcol(0.0001 .* q_Spin_Up_R1))) + 
    #     sum(sum(eachcol(0.0001 .* q_Spin_Up_R2))) + 
    #     sum(sum(eachcol(0.0001 .* q_Spin_Up_R3)))

    # reserve_payment = 
    #     sum(sum(eachcol(p_Spin_Up_R1 .* q_Spin_Up_R1)))

    # reserve_cost = 
    #     sum(sum(eachcol(0.0001 .* q_Spin_Up_R1)))

    # writedlm(joinpath(result_file_path, "reserve_payment.csv"), reserve_payment)
    # writedlm(joinpath(result_file_path, "reserve_cost.csv"), reserve_cost)

    # reserve_revenue_by_gen = DataFrame()
    # for g in get_components(Generator, sys)
    #     reserve_revenue_by_gen[!, get_name(g)] = zeros(time_step)
    # end

    # reserve_cost_by_gen = DataFrame()
    # for g in get_components(Generator, sys)
    #     reserve_cost_by_gen[!, get_name(g)] = zeros(time_step)
    # end

    # for g in get_components(Generator, sys)
    #     # if get_name(g) in names(q_Reg_Up)
    #     #     reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Reg_Up .* q_Reg_Up[!, get_name(g)]
    #     #     reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Reg_Up[!, get_name(g)]
    #     # end
    #     # if get_name(g) in names(q_Reg_Down)
    #     #     reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Reg_Down .* q_Reg_Down[!, get_name(g)]
    #     #     reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Reg_Down[!, get_name(g)]
    #     # end
    #     # if get_name(g) in names(q_Flex_Up)
    #     #     reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Flex_Up .* q_Flex_Up[!, get_name(g)]
    #     #     reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Flex_Up[!, get_name(g)]
    #     # end
    #     # if get_name(g) in names(q_Flex_Down)
    #     #     reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Flex_Down .* q_Flex_Down[!, get_name(g)]
    #     #     reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Flex_Down[!, get_name(g)]
    #     # end
    #     if get_name(g) in names(q_Spin_Up_R1)
    #         reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Spin_Up_R1 .* q_Spin_Up_R1[!, get_name(g)]
    #         reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Spin_Up_R1[!, get_name(g)]
    #     end
    #     # if get_name(g) in names(q_Spin_Up_R2)
    #     #     reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Spin_Up_R2 .* q_Spin_Up_R2[!, get_name(g)]
    #     #     reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Spin_Up_R2[!, get_name(g)]
    #     # end
    #     # if get_name(g) in names(q_Spin_Up_R3)
    #     #     reserve_revenue_by_gen[!, get_name(g)] = reserve_revenue_by_gen[!, get_name(g)] + p_Spin_Up_R3 .* q_Spin_Up_R3[!, get_name(g)]
    #     #     reserve_cost_by_gen[!, get_name(g)] = reserve_cost_by_gen[!, get_name(g)] + 0.0001 .* q_Spin_Up_R3[!, get_name(g)]
    #     # end
    # end

    ################### total load payment ########################
    total_load_payment = 0.0

    if load_output_model == "bid_in_load"
        load_payment_datacenter = DataFrame()
        for l in get_components(FlexibleLoad, sys)
            load_payment_datacenter[!, Symbol(get_name(l))] = zeros(time_step)
        end
        for l in get_components(FixedFlexibleLoad, sys, x -> PSY.get_available(x) == true)
            fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
            total_load_payment = total_load_payment + sum(LMP[1:time_step, get_name(l)] .* fixed_load)
        end
        for l in get_components(FlexibleLoad, sys)
            total_load_payment = total_load_payment + sum(LMP[!, get_name(l)] .* flex_load[!, get_name(l)])
        end
        for l in get_components(FlexibleLoad, sys)
            load_payment_datacenter[!, Symbol(get_name(l))] = LMP[!, get_name(l)] .* flex_load[!, get_name(l)]
        end
    elseif load_output_model == "constant_load"
        load_payment_datacenter = DataFrame()
        for l in get_components(PowerLoad, sys)
            load_payment_datacenter[!, Symbol(get_name(l))] = zeros(time_step)
        end
        for l in get_components(PowerLoad, sys)
            fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
            total_load_payment = total_load_payment + sum(LMP[1:time_step, get_name(l)] .* fixed_load)
            load_payment_datacenter[!, Symbol(get_name(l))] = LMP[!, get_name(l)] .* fixed_load
        end
    elseif load_output_model == "shiftable_load"
        load_payment_datacenter = DataFrame()
        for l in get_components(PowerLoad, sys)
            load_payment_datacenter[!, Symbol(get_name(l))] = zeros(time_step)
        end
        for l in get_components(PowerLoad, sys, x -> PSY.get_available(x) == true)
            fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
            total_load_payment = total_load_payment + sum(LMP[1:time_step, get_name(l)] .* fixed_load)
            load_payment_datacenter[!, Symbol(get_name(l))] = LMP[!, get_name(l)] .* fixed_load
        end
        for l in get_components(ShiftableLoad, sys)
            total_load_payment = total_load_payment + sum(LMP[1:time_step, get_name(l)] .* shift_load[1:time_step, get_name(l)])
        end
        for l in get_components(ShiftableLoad, sys)
            load_payment_datacenter[!, Symbol(get_name(l))] = load_payment_datacenter[!, Symbol(get_name(l))] .+ LMP[!, get_name(l)] .* shift_load[!, get_name(l)]
        end
    else
        error("Unexpected value: $load_output_model")
    end
    
    writedlm(joinpath(result_file_path, "total_load_payment.csv"), total_load_payment)
    CSV.write(joinpath(result_file_path, "load_payment_datacenter.csv"), load_payment_datacenter)

    ################### total load ########################
    total_load = 0.0

    if load_output_model == "bid_in_load"
        for l in get_components(FixedFlexibleLoad, sys, x -> PSY.get_available(x) == true)
            fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
            total_load = total_load + sum(fixed_load)
        end
        for l in get_components(FlexibleLoad, sys)
            total_load = total_load + sum(flex_load[!, get_name(l)])
        end
    elseif load_output_model == "constant_load"
        for l in get_components(PowerLoad, sys)
            fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
            total_load = total_load + sum(fixed_load)
        end
    elseif load_output_model == "shiftable_load"
        for l in get_components(PowerLoad, sys, x -> PSY.get_available(x) == true)
            fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
            total_load = total_load + sum(fixed_load)
        end
        for l in get_components(ShiftableLoad, sys)
            total_load = total_load + sum(shift_load[1:time_step, get_name(l)])
        end
    else
        error("Unexpected value: $load_output_model")
    end
    
    writedlm(joinpath(result_file_path, "total_load.csv"), total_load)

    load_fixed_ts = DataFrame()
    for l in get_components(FixedFlexibleLoad, sys)
        load_fixed_ts[!, Symbol(get_name(l))] = zeros(time_step)
    end
    for l in get_components(FixedFlexibleLoad, sys)
        fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[1:time_step] .* base_power
        load_fixed_ts[!, Symbol(get_name(l))] = fixed_load
    end

    CSV.write(joinpath(result_file_path, "load_fixed_ts.csv"), load_fixed_ts)

    ################### total value of load ########################
    if load_output_model == "bid_in_load"
        value_of_load = DataFrame()
        for l in get_components(FlexibleLoad, sys)
            value_of_load[!, Symbol(get_name(l))] = Float64[]
        end
        value_of_load_dict = Dict()
        for l in get_components(FlexibleLoad, sys)
            push!(value_of_load_dict, Symbol(get_name(l)) => 0.0)
        end
        push!(value_of_load, value_of_load_dict)

        for l in get_components(FlexibleLoad, sys)
            for t in 1:time_step
                flex_load_ts = variable_values["ActivePowerVariable__FlexibleLoad"][!, names(variable_values["ActivePowerVariable__FlexibleLoad"]) .== get_name(l)][t, 1]
                max_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                if flex_load_ts <= 1/5 * max_load
                    value_of_load[1, Symbol(get_name(l))] = value_of_load[1, Symbol(get_name(l))] + flex_load_ts * 1000
                elseif (flex_load_ts > 1/5 * max_load) & (flex_load_ts <= 2/5 * max_load)
                    value_of_load[1, Symbol(get_name(l))] = value_of_load[1, Symbol(get_name(l))] + max_load/5 * 1000 + (flex_load_ts-max_load/5) * 500
                elseif (flex_load_ts > 2/5 * max_load) & (flex_load_ts <= 3/5 * max_load)
                    value_of_load[1, Symbol(get_name(l))] = value_of_load[1, Symbol(get_name(l))] + max_load/5 * 1000 + max_load/5 * 500 + (flex_load_ts-max_load/5*2) * 100
                elseif (flex_load_ts > 3/5 * max_load) & (flex_load_ts <= 4/5 * max_load)
                    value_of_load[1, Symbol(get_name(l))] = value_of_load[1, Symbol(get_name(l))] + max_load/5 * 1000 + max_load/5 * 500 + max_load/5 * 100 + (flex_load_ts-max_load/5*3) * 50
                elseif flex_load_ts > 4/5 * max_load
                    value_of_load[1, Symbol(get_name(l))] = value_of_load[1, Symbol(get_name(l))] + max_load/5 * 1000 + max_load/5 * 500 + max_load/5 * 100 + max_load/5 * 50 + (flex_load_ts-max_load/5*4) * 10
                end
            end
        end
        CSV.write(joinpath(result_file_path, "value_of_load_by_data_center.csv"), value_of_load)
    end

    ###################### cost, value of load by time #########################
    if load_output_model == "bid_in_load"
        value_of_load_ts = DataFrame(load_value = zeros(time_step))
        for t in 1:time_step
            for l in get_components(FlexibleLoad, sys)
                flex_load_ts = variable_values["ActivePowerVariable__FlexibleLoad"][!, names(variable_values["ActivePowerVariable__FlexibleLoad"]) .== get_name(l)][t, 1]
                max_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                if flex_load_ts <= 1/5 * max_load
                    value_of_load_ts[t, :load_value] = value_of_load_ts[t, :load_value] + flex_load_ts * 1000
                elseif (flex_load_ts > 1/5 * max_load) & (flex_load_ts <= 2/5 * max_load)
                    value_of_load_ts[t, :load_value] = value_of_load_ts[t, :load_value] + max_load/5 * 1000 + (flex_load_ts-max_load/5) * 500
                elseif (flex_load_ts > 2/5 * max_load) & (flex_load_ts <= 3/5 * max_load)
                    value_of_load_ts[t, :load_value] = value_of_load_ts[t, :load_value] + max_load/5 * 1000 + max_load/5 * 500 + (flex_load_ts-max_load/5*2) * 100
                elseif (flex_load_ts > 3/5 * max_load) & (flex_load_ts <= 4/5 * max_load)
                    value_of_load_ts[t, :load_value] = value_of_load_ts[t, :load_value] + max_load/5 * 1000 + max_load/5 * 500 + max_load/5 * 100 + (flex_load_ts-max_load/5*3) * 50
                elseif flex_load_ts > 4/5 * max_load
                    value_of_load_ts[t, :load_value] = value_of_load_ts[t, :load_value] + max_load/5 * 1000 + max_load/5 * 500 + max_load/5 * 100 + max_load/5 * 50 + (flex_load_ts-max_load/5*4) * 10
                end
            end

            for l in get_components(FixedFlexibleLoad, sys, x -> PSY.get_available(x) == true)
                fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] * base_power
                value_of_load_ts[t, :load_value] = value_of_load_ts[t, :load_value] + 332 * fixed_load
            end
        end

        thermal_cost_ts = DataFrame(thermal_cost = zeros(time_step))
        for t in 1:time_step
            for l in get_components(ThermalStandard, sys)
                thermal_gen = variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(l)][t, 1]
                thermal_start = variable_values["StartVariable__ThermalStandard"][!, names(variable_values["StartVariable__ThermalStandard"]) .== get_name(l)][t, 1]
                thermal_shut = variable_values["StopVariable__ThermalStandard"][!, names(variable_values["StopVariable__ThermalStandard"]) .== get_name(l)][t, 1]
                thermal_on = variable_values["OnVariable__ThermalStandard"][!, names(variable_values["OnVariable__ThermalStandard"]) .== get_name(l)][t, 1]

                cost_start = thermal_start * get_start_up(get_operation_cost(l))
                cost_shut = thermal_shut * get_shut_down(get_operation_cost(l))
                cost_on = thermal_on * get_fixed(get_operation_cost(l))

                variable_cost_go = []
                variable_cost_psy = get_variable(get_operation_cost(l)).cost
                for b in 1:length(variable_cost_psy)
                    if b == 1
                        if length(variable_cost_psy[b]) == 1
                            # for this sys, this is the case of the "<1|2|3>14_SYNC_COND_1" devices with zero cost
                            block = [get_variable(get_operation_cost(l)).cost, get_max_active_power(l) * base_power]
                        else
                            block = [variable_cost_psy[b][1]/variable_cost_psy[b][2], variable_cost_psy[b][2]]
                        end
                    else
                        block = [(variable_cost_psy[b][1]-variable_cost_psy[b-1][1])/(variable_cost_psy[b][2]-variable_cost_psy[b-1][2]), (variable_cost_psy[b][2]-variable_cost_psy[b-1][2])]
                    end
                    push!(variable_cost_go, block)
                end

                cost_var = 0.0

                if length(variable_cost_psy) >= 2
                    # devices with # length(variable_cost_psy) ==1 are the "<1|2|3>14_SYNC_COND_1" devices with zero cost
                    for b in 1:length(variable_cost_psy)
                        if (b == 1) & (thermal_gen <= variable_cost_psy[b][2] + 0.000001)
                            cost_var = thermal_gen * variable_cost_go[b][1]
                        elseif b >= 2
                            if (thermal_gen > variable_cost_psy[b-1][2]) & (thermal_gen <= variable_cost_psy[b][2] + 0.000001)
                                cost_var = variable_cost_psy[b-1][1] + (thermal_gen-variable_cost_psy[b-1][2]) * variable_cost_go[b][1]
                            end
                        end
                    end
                end

                thermal_cost_ts[t, :thermal_cost] = thermal_cost_ts[t, :thermal_cost] + cost_start + cost_shut + cost_on + cost_var
            end
        end

        surplus_ts = DataFrame(surplus = value_of_load_ts[!, :load_value] - thermal_cost_ts[!, :thermal_cost])
        summary = DataFrame(surplus=surplus_ts[!, :surplus], load_value=value_of_load_ts[!, :load_value], cost=thermal_cost_ts[!, :thermal_cost])
        CSV.write(joinpath(result_file_path, "surplus_result_flex.csv"), summary)
    end

    ###################### revenue, cost, uplift and congestion rents #########################
    # generator revenue
    total_vg_revenue = 0.0
    total_non_vg_revenue = 0.0
    gen_revenue_total_ts = DataFrame()
    for g in get_components(Generator, sys)
        gen_revenue_total_ts[!, get_name(g)] = zeros(time_step)
    end

    for t in 1:time_step
        for g in get_components(Generator, sys)
            if typeof(g) == PSY.ThermalStandard
                dispatch = variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(g)][t, 1]
                total_non_vg_revenue = total_non_vg_revenue + LMP[t, get_name(get_bus(g))] * dispatch
            elseif typeof(g) == PSY.RenewableDispatch
                dispatch = variable_values["ActivePowerVariable__RenewableDispatch"][!, names(variable_values["ActivePowerVariable__RenewableDispatch"]) .== get_name(g)][t, 1]
                total_vg_revenue = total_vg_revenue + LMP[t, get_name(get_bus(g))] * dispatch
            elseif typeof(g) == PSY.HydroDispatch
                dispatch = variable_values["ActivePowerVariable__HydroDispatch"][!, names(variable_values["ActivePowerVariable__HydroDispatch"]) .== get_name(g)][t, 1]
                total_non_vg_revenue = total_non_vg_revenue + LMP[t, get_name(get_bus(g))] * dispatch
            elseif typeof(g) == PSY.HydroEnergyReservoir
                dispatch = variable_values["ActivePowerVariable__HydroEnergyReservoir"][!, names(variable_values["ActivePowerVariable__HydroEnergyReservoir"]) .== get_name(g)][t, 1]
                total_non_vg_revenue = total_non_vg_revenue + LMP[t, get_name(get_bus(g))] * dispatch
            elseif typeof(g) == PSY.RenewableFix
                ts = get_time_series_array(
                    SingleTimeSeries,
                    g,
                    "max_active_power",
                    start_time = start_date,
                    len = time_step
                )
                dispatch = TimeSeries.values(ts)[t] * base_power
                total_vg_revenue = total_vg_revenue + LMP[t, get_name(get_bus(g))] * dispatch
            else 
                error("Unexpected generator types")
            end

            gen_revenue_total_ts[t, get_name(g)] = LMP[t, get_name(get_bus(g))] * dispatch

        end
    end

    CSV.write(joinpath(result_file_path, "gen_revenue.csv"), gen_revenue_total_ts)
    writedlm(joinpath(result_file_path, "total_vg_revenue.csv"), total_vg_revenue)
    writedlm(joinpath(result_file_path, "total_non_vg_revenue.csv"), total_non_vg_revenue)

    # generator costs
    thermal_cost_ts = DataFrame(thermal_cost_variable = zeros(time_step), thermal_cost_fixed = zeros(time_step), thermal_cost_total = zeros(time_step))
    gen_cost_total_ts = DataFrame()
    for g in get_components(Generator, sys)
        gen_cost_total_ts[!, Symbol(get_name(g))] = zeros(time_step)
    end
    for t in 1:time_step
        for l in get_components(ThermalStandard, sys)
            thermal_gen = variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(l)][t, 1]
            thermal_start = variable_values["StartVariable__ThermalStandard"][!, names(variable_values["StartVariable__ThermalStandard"]) .== get_name(l)][t, 1]
            thermal_shut = variable_values["StopVariable__ThermalStandard"][!, names(variable_values["StopVariable__ThermalStandard"]) .== get_name(l)][t, 1]
            thermal_on = variable_values["OnVariable__ThermalStandard"][!, names(variable_values["OnVariable__ThermalStandard"]) .== get_name(l)][t, 1]

            cost_start = thermal_start * get_start_up(get_operation_cost(l))
            cost_shut = thermal_shut * get_shut_down(get_operation_cost(l))
            cost_on = thermal_on * get_fixed(get_operation_cost(l))

            variable_cost_go = []
            variable_cost_psy = get_variable(get_operation_cost(l)).cost
            for b in 1:length(variable_cost_psy)
                if b == 1
                    if length(variable_cost_psy[b]) == 1
                        # for this sys, this is the case of the "<1|2|3>14_SYNC_COND_1" devices with zero cost
                        block = [get_variable(get_operation_cost(l)).cost, get_max_active_power(l) * base_power]
                    else
                        block = [variable_cost_psy[b][1]/variable_cost_psy[b][2], variable_cost_psy[b][2]]
                    end
                else
                    block = [(variable_cost_psy[b][1]-variable_cost_psy[b-1][1])/(variable_cost_psy[b][2]-variable_cost_psy[b-1][2]), (variable_cost_psy[b][2]-variable_cost_psy[b-1][2])]
                end
                push!(variable_cost_go, block)
            end

            cost_var = 0.0

            if length(variable_cost_psy) >= 2
                # devices with # length(variable_cost_psy) ==1 are the "<1|2|3>14_SYNC_COND_1" devices with zero cost
                for b in 1:length(variable_cost_psy)
                    if (b == 1) & (thermal_gen <= variable_cost_psy[b][2] + 0.000001)
                        cost_var = thermal_gen * variable_cost_go[b][1]
                    elseif b >= 2
                        if (thermal_gen > variable_cost_psy[b-1][2]) & (thermal_gen <= variable_cost_psy[b][2] + 0.000001)
                            cost_var = variable_cost_psy[b-1][1] + (thermal_gen-variable_cost_psy[b-1][2]) * variable_cost_go[b][1]
                        end
                    end
                end
            end

            thermal_cost_ts[t, :thermal_cost_variable] = thermal_cost_ts[t, :thermal_cost_variable] + cost_var
            thermal_cost_ts[t, :thermal_cost_fixed] = thermal_cost_ts[t, :thermal_cost_fixed] + cost_start + cost_shut + cost_on
            thermal_cost_ts[t, :thermal_cost_total] = thermal_cost_ts[t, :thermal_cost_total] + cost_start + cost_shut + cost_on + cost_var
            gen_cost_total_ts[t, Symbol(get_name(l))] = cost_start + cost_shut + cost_on + cost_var
        end
    end

    CSV.write(joinpath(result_file_path, "thermal_cost.csv"), thermal_cost_ts)

    # uplift
    days = Int(time_step / 24)
    uplift_ts = DataFrame()
    for g in get_components(Generator, sys)
        uplift_ts[!, get_name(g)] = zeros(days)
    end
    uplift_vg = 0.0
    uplift_non_vg = 0.0

    for d in 1:days
        for g in get_components(Generator, sys)
            # uplift_ts[d, get_name(g)] = max(0.0, sum(gen_cost_total_ts[24*(d-1)+1:24*d, get_name(g)]) - sum(gen_revenue_total_ts[24*(d-1)+1:24*d, get_name(g)]) + 
            # sum(reserve_revenue_by_gen[24*(d-1)+1:24*d, get_name(g)]) - sum(reserve_cost_by_gen[24*(d-1)+1:24*d, get_name(g)]))
            uplift_ts[d, get_name(g)] = max(0.0, sum(gen_cost_total_ts[24*(d-1)+1:24*d, get_name(g)]) - sum(gen_revenue_total_ts[24*(d-1)+1:24*d, get_name(g)]))

            if typeof(g) == PSY.RenewableDispatch || typeof(g) == PSY.RenewableFix
                uplift_vg = uplift_vg + uplift_ts[d, get_name(g)]
            else
                uplift_non_vg = uplift_non_vg + uplift_ts[d, get_name(g)]
            end
        end
    end

    CSV.write(joinpath(result_file_path, "uplift.csv"), uplift_ts)
    writedlm(joinpath(result_file_path, "uplift_vg.csv"), uplift_vg)
    writedlm(joinpath(result_file_path, "uplift_non_vg.csv"), uplift_non_vg)

    # congestion rents
    load_payment_ts = DataFrame(load_payment = zeros(time_step))

    for t in 1:time_step
        if load_output_model == "bid_in_load"
            for l in get_components(FixedFlexibleLoad, sys, x -> PSY.get_available(x) == true)
                fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] * base_power
                load_payment_ts[t, :load_payment] = load_payment_ts[t, :load_payment] + LMP[t, get_name(l)] * fixed_load
            end
            for l in get_components(FlexibleLoad, sys)
                load_payment_ts[t, :load_payment] = load_payment_ts[t, :load_payment] + LMP[t, get_name(l)] * flex_load[t, get_name(l)]
            end
        elseif load_output_model == "constant_load"
            for l in get_components(PowerLoad, sys, x -> PSY.get_available(x) == true)
                fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] * base_power
                load_payment_ts[t, :load_payment] = load_payment_ts[t, :load_payment] + LMP[t, get_name(l)] * fixed_load
            end
        elseif load_output_model == "shiftable_load"
            for l in get_components(PowerLoad, sys, x -> PSY.get_available(x) == true)
                fixed_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] * base_power
                load_payment_ts[t, :load_payment] = load_payment_ts[t, :load_payment] + LMP[t, get_name(l)] * fixed_load
            end
            for l in get_components(ShiftableLoad, sys)
                load_payment_ts[t, :load_payment] = load_payment_ts[t, :load_payment] + LMP[t, get_name(l)] * shift_load[t, get_name(l)]
            end
        else
            error("Unexpected value: $load_output_model")
        end
    end

    congestion_rents_ts = DataFrame(congestion_rents = zeros(time_step))
    for t in 1:time_step
        congestion_rents_ts[t, :congestion_rents] = load_payment_ts[t, :load_payment] - sum(gen_revenue_total_ts[t, :])
    end

    CSV.write(joinpath(result_file_path, "congestion_rents.csv"), congestion_rents_ts)

    ###################### generation by fuel type #########################
    generation_by_fuel_ts = DataFrame(natural_gas = zeros(time_step), coal = zeros(time_step), oil = zeros(time_step), nuclear = zeros(time_step), other = zeros(time_step), solar = zeros(time_step), wind = zeros(time_step), hydro = zeros(time_step))
    for t in 1:time_step
        for g in get_components(ThermalStandard, sys)
            if string(get_fuel(g)) == "NATURAL_GAS"
                generation_by_fuel_ts[t, :natural_gas] = generation_by_fuel_ts[t, :natural_gas] + variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(g)][t, 1]
            elseif string(get_fuel(g)) == "COAL"
                generation_by_fuel_ts[t, :coal] = generation_by_fuel_ts[t, :coal] + variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(g)][t, 1]
            elseif string(get_fuel(g)) == "DISTILLATE_FUEL_OIL"
                generation_by_fuel_ts[t, :oil] = generation_by_fuel_ts[t, :oil] + variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(g)][t, 1]
            elseif string(get_fuel(g)) == "NUCLEAR"
                generation_by_fuel_ts[t, :nuclear] = generation_by_fuel_ts[t, :nuclear] + variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(g)][t, 1]
            else
                generation_by_fuel_ts[t, :other] = generation_by_fuel_ts[t, :other] + variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(g)][t, 1]
            end
        end
    
        for g in get_components(HydroDispatch, sys)
            generation_by_fuel_ts[t, :hydro] = generation_by_fuel_ts[t, :hydro] + variable_values["ActivePowerVariable__HydroDispatch"][!, names(variable_values["ActivePowerVariable__HydroDispatch"]) .== get_name(g)][t, 1]
        end
        for g in get_components(HydroEnergyReservoir, sys)
            generation_by_fuel_ts[t, :hydro] = generation_by_fuel_ts[t, :hydro] + variable_values["ActivePowerVariable__HydroEnergyReservoir"][!, names(variable_values["ActivePowerVariable__HydroEnergyReservoir"]) .== get_name(g)][t, 1]
        end
    
        for g in get_components(RenewableDispatch, sys)
            if get_prime_mover(g).value == 22
                generation_by_fuel_ts[t, :wind] = generation_by_fuel_ts[t, :wind] + variable_values["ActivePowerVariable__RenewableDispatch"][!, names(variable_values["ActivePowerVariable__RenewableDispatch"]) .== get_name(g)][t, 1]
            else
                generation_by_fuel_ts[t, :solar] = generation_by_fuel_ts[t, :solar] + variable_values["ActivePowerVariable__RenewableDispatch"][!, names(variable_values["ActivePowerVariable__RenewableDispatch"]) .== get_name(g)][t, 1]
            end
        end
    
        for g in get_components(RenewableFix, sys)
    
            ts = get_time_series_array(
                SingleTimeSeries,
                g,
                "max_active_power",
                start_time = start_date,
                len = time_step
            )
            ta = TimeSeries.values(ts) .* base_power
    
            if get_prime_mover(g).value == 22
                generation_by_fuel_ts[t, :wind] = generation_by_fuel_ts[t, :wind] + ta[t]
            else
                generation_by_fuel_ts[t, :solar] = generation_by_fuel_ts[t, :solar] + ta[t]
            end
        end
    
    end
    CSV.write(joinpath(result_file_path, "generation_by_fuel_ts.csv"), generation_by_fuel_ts)

    ###################### emission reporting #########################
    emission_total_ts = DataFrame()
    for e in emission_type
        emission_total_ts[!, Symbol(e)] = zeros(time_step)
    end
    for t in 1:time_step
        for l in get_components(ThermalStandard, sys)
            for e in emission_type
                thermal_gen = variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(l)][t, 1]
                thermal_start = variable_values["StartVariable__ThermalStandard"][!, names(variable_values["StartVariable__ThermalStandard"]) .== get_name(l)][t, 1]

                if source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1] != "Unit-specific"
                    if typeof(source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1]) == String || typeof(source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1]) == String15
                        emission_rate = parse(Float64, source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1])
                    else
                        emission_rate = source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1]
                    end
                else
                    emission_rate = 0.0
                end

                emission_start = thermal_start * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Start Heat Cold MBTU"][1] * emission_rate

                variable_cost_go = []
                variable_cost_psy = get_variable(get_operation_cost(l)).cost
                for b in 1:length(variable_cost_psy)
                    if b == 1
                        if length(variable_cost_psy[b]) == 1
                            block = [get_variable(get_operation_cost(l)).cost, get_max_active_power(l) * base_power]
                        else
                            block = [variable_cost_psy[b][1]/variable_cost_psy[b][2], variable_cost_psy[b][2]]
                        end
                    else
                        block = [(variable_cost_psy[b][1]-variable_cost_psy[b-1][1])/(variable_cost_psy[b][2]-variable_cost_psy[b-1][2]), (variable_cost_psy[b][2]-variable_cost_psy[b-1][2])]
                    end
                    push!(variable_cost_go, block)
                end

                emission_var = 0.0
                if length(variable_cost_psy) >= 2
                    if thermal_gen <= variable_cost_psy[1][2] + 0.000001
                        emission_var = thermal_gen * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate
                    elseif (thermal_gen > variable_cost_psy[1][2]) & (thermal_gen <= variable_cost_psy[2][2] + 0.000001)
                        emission_var = variable_cost_psy[1][2] * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate + 
                        (thermal_gen - variable_cost_psy[1][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_1"][1] / 1000 * emission_rate
                    elseif (thermal_gen > variable_cost_psy[2][2]) & (thermal_gen <= variable_cost_psy[3][2] + 0.000001)
                        emission_var = variable_cost_psy[1][2] * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate + 
                        (variable_cost_psy[2][2] - variable_cost_psy[1][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_1"][1] / 1000 * emission_rate + 
                        (thermal_gen - variable_cost_psy[2][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_2"][1] / 1000 * emission_rate
                    elseif (thermal_gen > variable_cost_psy[3][2]) & (thermal_gen <= variable_cost_psy[4][2] + 0.000001)
                        emission_var = variable_cost_psy[1][2] * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate + 
                        (variable_cost_psy[2][2] - variable_cost_psy[1][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_1"][1] / 1000 * emission_rate + 
                        (variable_cost_psy[3][2] - variable_cost_psy[2][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_2"][1] / 1000 * emission_rate + 
                        (thermal_gen - variable_cost_psy[3][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_3"][1] / 1000 * emission_rate
                    end
                end

                emission_total_ts[t, Symbol(e)] = emission_total_ts[t, Symbol(e)] + emission_start + emission_var
            end
        end
    end
    CSV.write(joinpath(result_file_path, "emission_total.csv"), emission_total_ts)

    ############# average emission rate #############
    total_emission_by_region = DataFrame()
    total_emission_by_region[!, Symbol("Region_1")] = zeros(time_step)
    total_emission_by_region[!, Symbol("Region_2")] = zeros(time_step)
    total_emission_by_region[!, Symbol("Region_3")] = zeros(time_step)
    for t in 1:time_step
        for l in get_components(ThermalStandard, sys)
            e = "CO2"
            thermal_gen = variable_values["ActivePowerVariable__ThermalStandard"][!, names(variable_values["ActivePowerVariable__ThermalStandard"]) .== get_name(l)][t, 1]
            thermal_start = variable_values["StartVariable__ThermalStandard"][!, names(variable_values["StartVariable__ThermalStandard"]) .== get_name(l)][t, 1]

            if source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1] != "Unit-specific"
                if typeof(source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1]) == String || typeof(source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1]) == String15
                    emission_rate = parse(Float64, source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1])
                else
                    emission_rate = source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Emissions " * e * " Lbs/MMBTU"][1]
                end
            else
                emission_rate = 0.0
            end

            emission_start = thermal_start * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "Start Heat Cold MBTU"][1] * emission_rate

            variable_cost_go = []
            variable_cost_psy = get_variable(get_operation_cost(l)).cost
            for b in 1:length(variable_cost_psy)
                if b == 1
                    if length(variable_cost_psy[b]) == 1
                        block = [get_variable(get_operation_cost(l)).cost, get_max_active_power(l) * base_power]
                    else
                        block = [variable_cost_psy[b][1]/variable_cost_psy[b][2], variable_cost_psy[b][2]]
                    end
                else
                    block = [(variable_cost_psy[b][1]-variable_cost_psy[b-1][1])/(variable_cost_psy[b][2]-variable_cost_psy[b-1][2]), (variable_cost_psy[b][2]-variable_cost_psy[b-1][2])]
                end
                push!(variable_cost_go, block)
            end

            emission_var = 0.0
            if length(variable_cost_psy) >= 2
                if thermal_gen <= variable_cost_psy[1][2] + 0.000001
                    emission_var = thermal_gen * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate
                elseif (thermal_gen > variable_cost_psy[1][2]) & (thermal_gen <= variable_cost_psy[2][2] + 0.000001)
                    emission_var = variable_cost_psy[1][2] * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate + 
                    (thermal_gen - variable_cost_psy[1][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_1"][1] / 1000 * emission_rate
                elseif (thermal_gen > variable_cost_psy[2][2]) & (thermal_gen <= variable_cost_psy[3][2] + 0.000001)
                    emission_var = variable_cost_psy[1][2] * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate + 
                    (variable_cost_psy[2][2] - variable_cost_psy[1][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_1"][1] / 1000 * emission_rate + 
                    (thermal_gen - variable_cost_psy[2][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_2"][1] / 1000 * emission_rate
                elseif (thermal_gen > variable_cost_psy[3][2]) & (thermal_gen <= variable_cost_psy[4][2] + 0.000001)
                    emission_var = variable_cost_psy[1][2] * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_avg_0"][1] / 1000 * emission_rate + 
                    (variable_cost_psy[2][2] - variable_cost_psy[1][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_1"][1] / 1000 * emission_rate + 
                    (variable_cost_psy[3][2] - variable_cost_psy[2][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_2"][1] / 1000 * emission_rate + 
                    (thermal_gen - variable_cost_psy[3][2]) * source_gen_data[source_gen_data[!, "GEN UID"] .== get_name(l), "HR_incr_3"][1] / 1000 * emission_rate
                end
            end

            total_emission_by_region[t, Symbol("Region_" * first(get_name(l), 1))] = total_emission_by_region[t, Symbol("Region_" * first(get_name(l), 1))] + emission_start + emission_var
        end
    end

    CSV.write(joinpath(result_file_path, "emission_total_by_region.csv"), total_emission_by_region)

    total_load_by_region = DataFrame()
    total_load_by_region[!, Symbol("Region_1")] = zeros(time_step)
    total_load_by_region[!, Symbol("Region_2")] = zeros(time_step)
    total_load_by_region[!, Symbol("Region_3")] = zeros(time_step)

    if load_output_model == "bid_in_load"
        for t in 1:time_step
            for l in get_components(FixedFlexibleLoad, sys, x -> PSY.get_available(x) == true)
                if get_bus(l).number < 200
                    total_load_by_region[t, Symbol("Region_1")] = total_load_by_region[t, Symbol("Region_1")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                elseif get_bus(l).number >= 300
                    total_load_by_region[t, Symbol("Region_3")] = total_load_by_region[t, Symbol("Region_3")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                else
                    total_load_by_region[t, Symbol("Region_2")] = total_load_by_region[t, Symbol("Region_2")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                end
            end
            for l in get_components(FlexibleLoad, sys)
                if get_bus(l).number < 200
                    total_load_by_region[t, Symbol("Region_1")] = total_load_by_region[t, Symbol("Region_1")] + variable_values["ActivePowerVariable__FlexibleLoad"][!, names(variable_values["ActivePowerVariable__FlexibleLoad"]) .== get_name(l)][t, 1]
                elseif get_bus(l).number >= 300
                    total_load_by_region[t, Symbol("Region_3")] = total_load_by_region[t, Symbol("Region_3")] + variable_values["ActivePowerVariable__FlexibleLoad"][!, names(variable_values["ActivePowerVariable__FlexibleLoad"]) .== get_name(l)][t, 1]
                else
                    total_load_by_region[t, Symbol("Region_2")] = total_load_by_region[t, Symbol("Region_2")] + variable_values["ActivePowerVariable__FlexibleLoad"][!, names(variable_values["ActivePowerVariable__FlexibleLoad"]) .== get_name(l)][t, 1]
                end
            end
        end
    elseif load_output_model == "constant_load"
        for t in 1:time_step
            for l in get_components(PowerLoad, sys)
                if get_bus(l).number < 200
                    total_load_by_region[t, Symbol("Region_1")] = total_load_by_region[t, Symbol("Region_1")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                elseif get_bus(l).number >= 300
                    total_load_by_region[t, Symbol("Region_3")] = total_load_by_region[t, Symbol("Region_3")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                else
                    total_load_by_region[t, Symbol("Region_2")] = total_load_by_region[t, Symbol("Region_2")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                end
            end
        end
    elseif load_output_model == "shiftable_load"
        for t in 1:time_step
            for l in get_components(PowerLoad, sys, x -> PSY.get_available(x) == true)
                if get_bus(l).number < 200
                    total_load_by_region[t, Symbol("Region_1")] = total_load_by_region[t, Symbol("Region_1")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                elseif get_bus(l).number >= 300
                    total_load_by_region[t, Symbol("Region_3")] = total_load_by_region[t, Symbol("Region_3")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                else
                    total_load_by_region[t, Symbol("Region_2")] = total_load_by_region[t, Symbol("Region_2")] + values(get_time_series_array(SingleTimeSeries, l, "max_active_power", start_time = start_date))[t] .* base_power
                end
            end
            for l in get_components(ShiftableLoad, sys)
                if get_bus(l).number < 200
                    total_load_by_region[t, Symbol("Region_1")] = total_load_by_region[t, Symbol("Region_1")] + shift_load[t, get_name(l)]
                elseif get_bus(l).number >= 300
                    total_load_by_region[t, Symbol("Region_3")] = total_load_by_region[t, Symbol("Region_3")] + shift_load[t, get_name(l)]
                else
                    total_load_by_region[t, Symbol("Region_2")] = total_load_by_region[t, Symbol("Region_2")] + shift_load[t, get_name(l)]
                end
            end
        end

    end

    CSV.write(joinpath(result_file_path, "load_total_by_region.csv"), total_load_by_region)

    average_emission_rate_by_region = total_emission_by_region ./ total_load_by_region

    CSV.write(joinpath(result_file_path, "average_emission_rate_by_region.csv"), average_emission_rate_by_region)

    return average_emission_rate_by_region, time_step
end