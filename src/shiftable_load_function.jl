const SECONDS_IN_HOUR = 3600.0

############# shiftable load formulation #############
mutable struct ShiftableLoad <: PSY.ControllableLoad
    name::String
    available::Bool
    bus::Bus
    model::PSY.LoadModels
    active_power::Float64
    reactive_power::Float64
    max_active_power::Float64
    max_reactive_power::Float64
    active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}}
    "Base power of the unit in MVA"
    base_power::Float64
    "Services that this device contributes to"
    services::Vector{PSY.Service}
    "corresponding dynamic injection device"
    dynamic_injector::Union{Nothing, PSY.DynamicInjection}
    ext::Dict{String, Any}
    "internal time_series storage"
    time_series_container::InfrastructureSystems.TimeSeriesContainer
    "power system internal reference, do not modify"
    internal::InfrastructureSystems.InfrastructureSystemsInternal
end

function ShiftableLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), )
    ShiftableLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, services, dynamic_injector, ext, time_series_container, InfrastructureSystemsInternal(), )
end

function ShiftableLoad(; name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), internal=InfrastructureSystems.InfrastructureSystemsInternal(), )
    ShiftableLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, services, dynamic_injector, ext, time_series_container, internal, )
end

"""Get [`InterruptibleLoad`](@ref) `name`."""
PSY.get_name(value::ShiftableLoad) = value.name
"""Get [`InterruptibleLoad`](@ref) `available`."""
PSY.get_available(value::ShiftableLoad) = value.available
"""Get [`InterruptibleLoad`](@ref) `bus`."""
PSY.get_bus(value::ShiftableLoad) = value.bus
"""Get [`InterruptibleLoad`](@ref) `model`."""
PSY.get_model(value::ShiftableLoad) = value.model
"""Get [`InterruptibleLoad`](@ref) `active_power`."""
PSY.get_active_power(value::ShiftableLoad) = PSY.get_value(value, value.active_power)
"""Get [`InterruptibleLoad`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::ShiftableLoad) = PSY.get_value(value, value.reactive_power)
"""Get [`InterruptibleLoad`](@ref) `max_active_power`."""
PSY.get_max_active_power(value::ShiftableLoad) = PSY.get_value(value, value.max_active_power)
"""Get [`InterruptibleLoad`](@ref) `max_reactive_power`."""
PSY.get_max_reactive_power(value::ShiftableLoad) = PSY.get_value(value, value.max_reactive_power)
"""Get [`InterruptibleLoad`](@ref) `base_power`."""
PSY.get_base_power(value::ShiftableLoad) = value.base_power
"""Get [`InterruptibleLoad`](@ref) `services`."""
PSY.get_services(value::ShiftableLoad) = value.services
"""Get [`InterruptibleLoad`](@ref) `dynamic_injector`."""
PSY.get_dynamic_injector(value::ShiftableLoad) = value.dynamic_injector
"""Get [`InterruptibleLoad`](@ref) `ext`."""
PSY.get_ext(value::ShiftableLoad) = value.ext
"""Get [`InterruptibleLoad`](@ref) `time_series_container`."""
PSY.get_time_series_container(value::ShiftableLoad) = value.time_series_container
"""Get [`InterruptibleLoad`](SECONDS_IN_HOUR@ref) `internal`."""
PSY.get_internal(value::ShiftableLoad) = value.internal
"""Get [`HydroDispatch`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::ShiftableLoad) = PSY.get_value(value, value.active_power_limits)

PSY.set_available!(value::ShiftableLoad, val) = value.available = val

struct LoadShifting <: PSI.AbstractControllablePowerLoadFormulation end

struct LoadBudgetConstraint <: PSI.ConstraintType end

function PSI.get_default_time_series_names(
    ::Type{<:ShiftableLoad},
    ::Type{<:LoadShifting},
)
    return Dict{Type{<:PSI.TimeSeriesParameter}, String}(
        PSI.ActivePowerTimeSeriesParameter => "max_active_power",
        PSI.EnergyBudgetTimeSeriesParameter => "load_budget",
    )
end


function PSI.add_constraints!(
    container::PSI.OptimizationContainer,
    T::Type{PSI.ActivePowerVariableLimitsConstraint},
    U::Type{<:Union{PSI.VariableType, PSI.ExpressionType}},
    devices::IS.FlattenIteratorWrapper{V},
    model::PSI.DeviceModel{V, W},
    X::Type{<:PM.AbstractPowerModel},
) where {V <: Union{ShiftableLoad, ChlorAlkaliPlantWithBidLoad}, W <: Union{LoadShifting, ChlorAlkaliPlantDispatchWithBid}}
    if !PSI.has_semicontinuous_feedforward(model, U)
        PSI.add_range_constraints!(container, T, U, devices, model, X)
    end
    PSI.add_parameterized_upper_bound_range_constraints(
        container,
        PSI.ActivePowerVariableTimeSeriesLimitsConstraint,
        U,
        PSI.ActivePowerTimeSeriesParameter,
        devices,
        model,
        X,
    )
    return
end


function PSI.add_constraints!(
    container::PSI.OptimizationContainer,
    ::Type{LoadBudgetConstraint},
    devices::IS.FlattenIteratorWrapper{V},
    model::PSI.DeviceModel{V, W},
    ::Type{X},
) where {V <: Union{ShiftableLoad, ChlorAlkaliPlantWithBidLoad}, W <: Union{LoadShifting, ChlorAlkaliPlantDispatchWithBid}, X <: PM.AbstractPowerModel}
    time_steps = PSI.get_time_steps(container)
    resolution = PSI.get_resolution(container)
    set_name = [PSY.get_name(d) for d in devices]
    constraint =
        PSI.add_constraints_container!(container, LoadBudgetConstraint(), V, set_name)

    variable_out = PSI.get_variable(container, PSI.ActivePowerVariable(), V)
    param = PSI.get_parameter_array(container, PSI.EnergyBudgetTimeSeriesParameter(), V)
    multiplier =
        PSI.get_parameter_multiplier_array(container, PSI.EnergyBudgetTimeSeriesParameter(), V)

    for d in devices
        name = PSY.get_name(d)
        constraint[name] = JuMP.@constraint(
            container.JuMPmodel,
            sum([variable_out[name, t] for t in time_steps]) == sum([multiplier[name, t] * param[name, t] for t in time_steps])
        )
    end
    return
end


function PSI.construct_device!(
    container::PSI.OptimizationContainer,
    sys::PSY.System,
    ::PSI.ArgumentConstructStage,
    model::PSI.DeviceModel{H, LoadShifting},
    ::Type{S},
) where {H <: ShiftableLoad, S <: PM.AbstractActivePowerModel}
    devices = PSI.get_available_components(H, sys)

    PSI.add_variables!(container, PSI.ActivePowerVariable, devices, LoadShifting())
    PSI.add_parameters!(container, PSI.ActivePowerTimeSeriesParameter, devices, model)
    PSI.add_parameters!(container, PSI.EnergyBudgetTimeSeriesParameter, devices, model)

    PSI.add_expressions!(container, PSI.ProductionCostExpression, devices, model)

    PSI.add_to_expression!(
        container,
        PSI.ActivePowerBalance,
        PSI.ActivePowerVariable,
        devices,
        model,
        S,
    )

    return
end

function PSI.construct_device!(
    container::PSI.OptimizationContainer,
    sys::PSY.System,
    ::PSI.ModelConstructStage,
    model::DeviceModel{H, LoadShifting},
    ::Type{S},
) where {H <: ShiftableLoad, S <: PM.AbstractActivePowerModel}
    devices = PSI.get_available_components(H, sys)

    PSI.add_constraints!(
        container,
        PSI.ActivePowerVariableLimitsConstraint,
        PSI.ActivePowerVariable,
        devices,
        model,
        S,
    )
    
    # Energy Budget Constraint
    PSI.add_constraints!(container, LoadBudgetConstraint, devices, model, S)

    PSI.add_feedforward_constraints!(container, model, devices)

    PSI.add_constraint_dual!(container, sys, model)

    return
end

function PSI.get_min_max_limits(
    x::ShiftableLoad,
    ::Type{<:PSI.ActivePowerVariableLimitsConstraint},
    ::Type{LoadShifting},
)
    return PSY.get_active_power_limits(x)
end






















# function PSI.DeviceRangeConstraintSpec(
#     ::Type{<:PSI.RangeConstraint},
#     ::Type{PSI.ActivePowerVariable},
#     ::Type{T},
#     ::Type{<:LoadShifting},
#     ::Type{<:PM.AbstractPowerModel},
#     feedforward::Union{Nothing, PSI.AbstractAffectFeedForward},
#     use_parameters::Bool,
#     use_forecasts::Bool,
# ) where {T <: PSY.ControllableLoad}
#     if !use_parameters && !use_forecasts
#         return PSI.DeviceRangeConstraintSpec(;
#             range_constraint_spec = PSI.RangeConstraintSpec(;
#                 constraint_name = PSI.make_constraint_name(
#                     PSI.RangeConstraint,
#                     PSI.ActivePowerVariable,
#                     T,
#                 ),
#                 variable_name = PSI.make_variable_name(PSI.ActivePowerVariable, T),
#                 limits_func = x -> (min = 0.0, max = PSY.get_active_power(x)),
#                 constraint_func = PSI.device_range!,
#                 constraint_struct = PSI.DeviceRangeConstraintInfo,
#             ),
#         )
#     end

#     return PSI.DeviceRangeConstraintSpec(;
#         timeseries_range_constraint_spec = PSI.TimeSeriesConstraintSpec(
#             constraint_name = PSI.make_constraint_name(
#                 PSI.RangeConstraint,
#                 PSI.ActivePowerVariable,
#                 T,
#             ),
#             variable_name = PSI.make_variable_name(PSI.ActivePowerVariable, T),
#             parameter_name = use_parameters ? PSI.ACTIVE_POWER : nothing,
#             forecast_label = "max_hourly_load_budget",
#             multiplier_func = x -> 1.0,
#             constraint_func = use_parameters ? PSI.device_timeseries_param_ub! :
#                               PSI.device_timeseries_ub!,
#         ),
#     )
# end

# function device_energy_budget_param_fb(
#     optimization_container::PSI.OptimizationContainer,
#     energy_budget_data::Vector{PSI.DeviceTimeSeriesConstraintInfo},
#     cons_name::Symbol,
#     param_reference::PSI.UpdateRef,
#     var_names::Symbol,
# )
#     time_steps = PSI.model_time_steps(optimization_container)
#     resolution = PSI.model_resolution(optimization_container)
#     inv_dt = 1.0 / (Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR)
#     variable_out = PSI.get_variable(optimization_container, var_names)
#     set_name = [PSI.get_component_name(r) for r in energy_budget_data]
#     constraint = PSI.add_cons_container!(optimization_container, cons_name, set_name)
#     container =
#         PSI.add_param_container!(optimization_container, param_reference, set_name, time_steps)
#     multiplier = PSI.get_multiplier_array(container)
#     param = PSI.get_parameter_array(container)
#     for constraint_info in energy_budget_data
#         name = PSI.get_component_name(constraint_info)
#         for t in time_steps
#             multiplier[name, t] = constraint_info.multiplier * inv_dt
#             param[name, t] = PSI.add_parameter(
#                 optimization_container.JuMPmodel,
#                 constraint_info.timeseries[t],
#             )
#         end
#         constraint[name] = JuMP.@constraint(
#             optimization_container.JuMPmodel,
#             sum([variable_out[name, t] for t in time_steps]) == sum([multiplier[name, t] * param[name, t] for t in time_steps])
#         )
#     end

#     return
# end

# function device_energy_budget_fb(
#     optimization_container::PSI.OptimizationContainer,
#     energy_budget_constraints::Vector{PSI.DeviceTimeSeriesConstraintInfo},
#     cons_name::Symbol,
#     var_names::Symbol,
# )
#     time_steps = PSI.model_time_steps(optimization_container)
#     variable_out = PSI.get_variable(optimization_container, var_names)
#     names = [PSI.get_component_name(x) for x in energy_budget_constraints]
#     constraint = PSI.add_cons_container!(optimization_container, cons_name, names)

#     for constraint_info in energy_budget_constraints
#         name = PSI.get_component_name(constraint_info)
#         resolution = PSI.model_resolution(optimization_container)
#         inv_dt = 1.0 / (Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR)
#         forecast = constraint_info.timeseries
#         multiplier = constraint_info.multiplier * inv_dt
#         constraint[name] = JuMP.@constraint(
#             optimization_container.JuMPmodel,
#             sum([variable_out[name, t] for t in time_steps]) == multiplier * sum(forecast)
#         )
#     end

#     return
# end

# function PSI.energy_budget_constraints!(
#     optimization_container::PSI.OptimizationContainer,
#     devices::IS.FlattenIteratorWrapper{H},
#     ::DeviceModel{H, <:LoadShifting},
#     ::Type{<:PM.AbstractPowerModel},
#     ::Union{Nothing, PSI.AbstractAffectFeedForward},
# ) where {H <: PSY.ControllableLoad}
#     forecast_label = "load_budget"
#     constraint_data = Vector{PSI.DeviceTimeSeriesConstraintInfo}(undef, length(devices))
#     for (ix, d) in enumerate(devices)
#         ts_vector = PSI.get_time_series(optimization_container, d, forecast_label)
#         @debug "time_series" ts_vector
#         # constraint_d =
#         #     PSI.DeviceTimeSeriesConstraintInfo(d, x -> PSY.get_storage_capacity(x), ts_vector)
#         constraint_d =
#             PSI.DeviceTimeSeriesConstraintInfo(d, x -> 1.0, ts_vector)
#         constraint_data[ix] = constraint_d
#     end

#     if PSI.model_has_parameters(optimization_container)
#         device_energy_budget_param_fb(
#             optimization_container,
#             constraint_data,
#             PSI.make_constraint_name(ENERGY_BUDGET, H),
#             PSI.UpdateRef{H}(ENERGY_BUDGET, forecast_label),
#             PSI.make_variable_name(ACTIVE_POWER, H),
#         )
#     else
#         device_energy_budget_fb(
#             optimization_container,
#             constraint_data,
#             PSI.make_constraint_name(ENERGY_BUDGET),
#             PSI.make_variable_name(ACTIVE_POWER, H),
#         )
#     end
# end

# function PSI.construct_device!(
#     optimization_container::PSI.OptimizationContainer,
#     sys::PSY.System,
#     model::PSI.DeviceModel{H, D},
#     ::Type{S},
# ) where {
#     H <: ShiftableLoad,
#     D <: LoadShifting,
#     S <: PM.AbstractActivePowerModel,
# }
#     devices = PSI.get_available_components(H, sys)

#     if !PSI.validate_available_devices(H, devices)
#         return
#     end

#     # Variables
#     PSI.add_variables!(optimization_container, PSI.ActivePowerVariable, devices, D())

#     # Energy Budget Constraint
#     PSI.energy_budget_constraints!(
#         optimization_container,
#         devices,
#         model,
#         S,
#         PSI.get_feedforward(model),
#     )

#     # Range Constraints
#     PSI.add_constraints!(
#         optimization_container,
#         PSI.RangeConstraint,
#         PSI.ActivePowerVariable,
#         devices,
#         model,
#         S,
#         PSI.get_feedforward(model),
#     )

#     PSI.feedforward!(optimization_container, devices, model, PSI.get_feedforward(model))

#     # Cost Function
#     # PSI.cost_function!(optimization_container, devices, model, S, nothing)

#     return
# end