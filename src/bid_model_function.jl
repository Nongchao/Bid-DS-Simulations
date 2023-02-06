const SECONDS_IN_HOUR = 3600.0

mutable struct FlexibleLoad <: PSY.ControllableLoad
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
    "Operation Cost of Generation [`TwoPartCost`](@ref)"
    operation_cost::PSY.OperationalCost
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

function FlexibleLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), )
    FlexibleLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, InfrastructureSystemsInternal(), )
end

function FlexibleLoad(; name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), internal=InfrastructureSystems.InfrastructureSystemsInternal(), )
    FlexibleLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, internal, )
end

"""Get [`InterruptibleLoad`](@ref) `name`."""
PSY.get_name(value::FlexibleLoad) = value.name
"""Get [`InterruptibleLoad`](@ref) `available`."""
PSY.get_available(value::FlexibleLoad) = value.available
"""Get [`InterruptibleLoad`](@ref) `bus`."""
PSY.get_bus(value::FlexibleLoad) = value.bus
"""Get [`InterruptibleLoad`](@ref) `model`."""
PSY.get_model(value::FlexibleLoad) = value.model
"""Get [`InterruptibleLoad`](@ref) `active_power`."""
PSY.get_active_power(value::FlexibleLoad) = PSY.get_value(value, value.active_power)
"""Get [`InterruptibleLoad`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::FlexibleLoad) = PSY.get_value(value, value.reactive_power)
"""Get [`InterruptibleLoad`](@ref) `max_active_power`."""
PSY.get_max_active_power(value::FlexibleLoad) = PSY.get_value(value, value.max_active_power)
"""Get [`InterruptibleLoad`](@ref) `max_reactive_power`."""
PSY.get_max_reactive_power(value::FlexibleLoad) = PSY.get_value(value, value.max_reactive_power)
"""Get [`InterruptibleLoad`](@ref) `base_power`."""
PSY.get_base_power(value::FlexibleLoad) = value.base_power
"""Get [`InterruptibleLoad`](@ref) `operation_cost`."""
PSY.get_operation_cost(value::FlexibleLoad) = value.operation_cost
"""Get [`InterruptibleLoad`](@ref) `services`."""
PSY.get_services(value::FlexibleLoad) = value.services
"""Get [`InterruptibleLoad`](@ref) `dynamic_injector`."""
PSY.get_dynamic_injector(value::FlexibleLoad) = value.dynamic_injector
"""Get [`InterruptibleLoad`](@ref) `ext`."""
PSY.get_ext(value::FlexibleLoad) = value.ext
"""Get [`InterruptibleLoad`](@ref) `time_series_container`."""
PSY.get_time_series_container(value::FlexibleLoad) = value.time_series_container
"""Get [`InterruptibleLoad`](@ref) `internal`."""
PSY.get_internal(value::FlexibleLoad) = value.internal
"""Get [`HydroDispatch`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::FlexibleLoad) = PSY.get_value(value, value.active_power_limits)

PSY.set_available!(value::FlexibleLoad, val) = value.available = val

mutable struct FixedFlexibleLoad <: PSY.ControllableLoad
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
    "Operation Cost of Generation [`TwoPartCost`](@ref)"
    operation_cost::PSY.OperationalCost
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

function FixedFlexibleLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), )
    FixedFlexibleLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, InfrastructureSystemsInternal(), )
end

function FixedFlexibleLoad(; name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), internal=InfrastructureSystems.InfrastructureSystemsInternal(), )
    FixedFlexibleLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, internal, )
end

"""Get [`InterruptibleLoad`](@ref) `name`."""
PSY.get_name(value::FixedFlexibleLoad) = value.name
"""Get [`InterruptibleLoad`](@ref) `available`."""
PSY.get_available(value::FixedFlexibleLoad) = value.available
"""Get [`InterruptibleLoad`](@ref) `bus`."""
PSY.get_bus(value::FixedFlexibleLoad) = value.bus
"""Get [`InterruptibleLoad`](@ref) `model`."""
PSY.get_model(value::FixedFlexibleLoad) = value.model
"""Get [`InterruptibleLoad`](@ref) `active_power`."""
PSY.get_active_power(value::FixedFlexibleLoad) = PSY.get_value(value, value.active_power)
"""Get [`InterruptibleLoad`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::FixedFlexibleLoad) = PSY.get_value(value, value.reactive_power)
"""Get [`InterruptibleLoad`](@ref) `max_active_power`."""
PSY.get_max_active_power(value::FixedFlexibleLoad) = PSY.get_value(value, value.max_active_power)
"""Get [`InterruptibleLoad`](@ref) `max_reactive_power`."""
PSY.get_max_reactive_power(value::FixedFlexibleLoad) = PSY.get_value(value, value.max_reactive_power)
"""Get [`InterruptibleLoad`](@ref) `base_power`."""
PSY.get_base_power(value::FixedFlexibleLoad) = value.base_power
"""Get [`InterruptibleLoad`](@ref) `operation_cost`."""
PSY.get_operation_cost(value::FixedFlexibleLoad) = value.operation_cost
"""Get [`InterruptibleLoad`](@ref) `services`."""
PSY.get_services(value::FixedFlexibleLoad) = value.services
"""Get [`InterruptibleLoad`](@ref) `dynamic_injector`."""
PSY.get_dynamic_injector(value::FixedFlexibleLoad) = value.dynamic_injector
"""Get [`InterruptibleLoad`](@ref) `ext`."""
PSY.get_ext(value::FixedFlexibleLoad) = value.ext
"""Get [`InterruptibleLoad`](@ref) `time_series_container`."""
PSY.get_time_series_container(value::FixedFlexibleLoad) = value.time_series_container
"""Get [`InterruptibleLoad`](@ref) `internal`."""
PSY.get_internal(value::FixedFlexibleLoad) = value.internal
"""Get [`HydroDispatch`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::FixedFlexibleLoad) = PSY.get_value(value, value.active_power_limits)

PSY.set_available!(value::FixedFlexibleLoad, val) = value.available = val

mutable struct CurtailableLoad <: PSY.ControllableLoad
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
    "Operation Cost of Generation [`TwoPartCost`](@ref)"
    operation_cost::PSY.OperationalCost
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

function CurtailableLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), )
    CurtailableLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, InfrastructureSystemsInternal(), )
end

function CurtailableLoad(; name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), internal=InfrastructureSystems.InfrastructureSystemsInternal(), )
    CurtailableLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, internal, )
end

"""Get [`InterruptibleLoad`](@ref) `name`."""
PSY.get_name(value::CurtailableLoad) = value.name
"""Get [`InterruptibleLoad`](@ref) `available`."""
PSY.get_available(value::CurtailableLoad) = value.available
"""Get [`InterruptibleLoad`](@ref) `bus`."""
PSY.get_bus(value::CurtailableLoad) = value.bus
"""Get [`InterruptibleLoad`](@ref) `model`."""
PSY.get_model(value::CurtailableLoad) = value.model
"""Get [`InterruptibleLoad`](@ref) `active_power`."""
PSY.get_active_power(value::CurtailableLoad) = PSY.get_value(value, value.active_power)
"""Get [`InterruptibleLoad`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::CurtailableLoad) = PSY.get_value(value, value.reactive_power)
"""Get [`InterruptibleLoad`](@ref) `max_active_power`."""
PSY.get_max_active_power(value::CurtailableLoad) = PSY.get_value(value, value.max_active_power)
"""Get [`InterruptibleLoad`](@ref) `max_reactive_power`."""
PSY.get_max_reactive_power(value::CurtailableLoad) = PSY.get_value(value, value.max_reactive_power)
"""Get [`InterruptibleLoad`](@ref) `base_power`."""
PSY.get_base_power(value::CurtailableLoad) = value.base_power
"""Get [`InterruptibleLoad`](@ref) `operation_cost`."""
PSY.get_operation_cost(value::CurtailableLoad) = value.operation_cost
"""Get [`InterruptibleLoad`](@ref) `services`."""
PSY.get_services(value::CurtailableLoad) = value.services
"""Get [`InterruptibleLoad`](@ref) `dynamic_injector`."""
PSY.get_dynamic_injector(value::CurtailableLoad) = value.dynamic_injector
"""Get [`InterruptibleLoad`](@ref) `ext`."""
PSY.get_ext(value::CurtailableLoad) = value.ext
"""Get [`InterruptibleLoad`](@ref) `time_series_container`."""
PSY.get_time_series_container(value::CurtailableLoad) = value.time_series_container
"""Get [`InterruptibleLoad`](@ref) `internal`."""
PSY.get_internal(value::CurtailableLoad) = value.internal
"""Get [`HydroDispatch`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::CurtailableLoad) = PSY.get_value(value, value.active_power_limits)

PSY.set_available!(value::CurtailableLoad, val) = value.available = val

mutable struct ChlorAlkaliPlantWithBidLoad <: PSY.ControllableLoad
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
    "Operation Cost of Generation [`TwoPartCost`](@ref)"
    operation_cost::PSY.OperationalCost
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

function ChlorAlkaliPlantWithBidLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), )
    ChlorAlkaliPlantWithBidLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, InfrastructureSystemsInternal(), )
end

function ChlorAlkaliPlantWithBidLoad(; name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services=Device[], dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), internal=InfrastructureSystems.InfrastructureSystemsInternal(), )
    ChlorAlkaliPlantWithBidLoad(name, available, bus, model, active_power, reactive_power, max_active_power, max_reactive_power, active_power_limits, base_power, operation_cost, services, dynamic_injector, ext, time_series_container, internal, )
end

"""Get [`InterruptibleLoad`](@ref) `name`."""
PSY.get_name(value::ChlorAlkaliPlantWithBidLoad) = value.name
"""Get [`InterruptibleLoad`](@ref) `available`."""
PSY.get_available(value::ChlorAlkaliPlantWithBidLoad) = value.available
"""Get [`InterruptibleLoad`](@ref) `bus`."""
PSY.get_bus(value::ChlorAlkaliPlantWithBidLoad) = value.bus
"""Get [`InterruptibleLoad`](@ref) `model`."""
PSY.get_model(value::ChlorAlkaliPlantWithBidLoad) = value.model
"""Get [`InterruptibleLoad`](@ref) `active_power`."""
PSY.get_active_power(value::ChlorAlkaliPlantWithBidLoad) = PSY.get_value(value, value.active_power)
"""Get [`InterruptibleLoad`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::ChlorAlkaliPlantWithBidLoad) = PSY.get_value(value, value.reactive_power)
"""Get [`InterruptibleLoad`](@ref) `max_active_power`."""
PSY.get_max_active_power(value::ChlorAlkaliPlantWithBidLoad) = PSY.get_value(value, value.max_active_power)
"""Get [`InterruptibleLoad`](@ref) `max_reactive_power`."""
PSY.get_max_reactive_power(value::ChlorAlkaliPlantWithBidLoad) = PSY.get_value(value, value.max_reactive_power)
"""Get [`InterruptibleLoad`](@ref) `base_power`."""
PSY.get_base_power(value::ChlorAlkaliPlantWithBidLoad) = value.base_power
"""Get [`InterruptibleLoad`](@ref) `operation_cost`."""
PSY.get_operation_cost(value::ChlorAlkaliPlantWithBidLoad) = value.operation_cost
"""Get [`InterruptibleLoad`](@ref) `services`."""
PSY.get_services(value::ChlorAlkaliPlantWithBidLoad) = value.services
"""Get [`InterruptibleLoad`](@ref) `dynamic_injector`."""
PSY.get_dynamic_injector(value::ChlorAlkaliPlantWithBidLoad) = value.dynamic_injector
"""Get [`InterruptibleLoad`](@ref) `ext`."""
PSY.get_ext(value::ChlorAlkaliPlantWithBidLoad) = value.ext
"""Get [`InterruptibleLoad`](@ref) `time_series_container`."""
PSY.get_time_series_container(value::ChlorAlkaliPlantWithBidLoad) = value.time_series_container
"""Get [`InterruptibleLoad`](@ref) `internal`."""
PSY.get_internal(value::ChlorAlkaliPlantWithBidLoad) = value.internal
"""Get [`HydroDispatch`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::ChlorAlkaliPlantWithBidLoad) = PSY.get_value(value, value.active_power_limits)

PSY.set_available!(value::ChlorAlkaliPlantWithBidLoad, val) = value.available = val


struct FixDispatchablePowerLoad <: PSI.AbstractControllablePowerLoadFormulation end
struct FixDispatchablePowerLoad_60perc <: PSI.AbstractControllablePowerLoadFormulation end
struct CurtailablePowerLoad <: PSI.AbstractControllablePowerLoadFormulation end
struct ChlorAlkaliPlantDispatchWithBid <: PSI.AbstractControllablePowerLoadFormulation end

PSI.objective_function_multiplier(::PSI.VariableType, ::FixDispatchablePowerLoad)=PSI.OBJECTIVE_FUNCTION_POSITIVE
PSI.objective_function_multiplier(::PSI.VariableType, ::FixDispatchablePowerLoad_60perc)=PSI.OBJECTIVE_FUNCTION_POSITIVE
PSI.objective_function_multiplier(::PSI.VariableType, ::PSI.DispatchablePowerLoad)=PSI.OBJECTIVE_FUNCTION_POSITIVE
PSI.objective_function_multiplier(::PSI.VariableType, ::CurtailablePowerLoad)=PSI.OBJECTIVE_FUNCTION_POSITIVE
PSI.objective_function_multiplier(::PSI.VariableType, ::ChlorAlkaliPlantDispatchWithBid)=PSI.OBJECTIVE_FUNCTION_POSITIVE

function fixed_bound_range_with_parameter!(
    jump_model::JuMP.Model,
    constraint_container::PSI.JuMPConstraintArray,
    lhs_array,
    param_multiplier::PSI.JuMPFloatArray,
    param_array::Union{PSI.JuMPParamArray, PSI.JuMPFloatArray},
    devices::IS.FlattenIteratorWrapper{V},
) where {V <: PSY.Component}
    time_steps = axes(constraint_container)[2]
    for device in devices, t in time_steps
        name = PSY.get_name(device)
        constraint_container[name, t] = JuMP.@constraint(
            jump_model,
            lhs_array[name, t] == param_multiplier[name, t] * param_array[name, t]
        )
    end
    return
end

function _add_parameterized_fixed_bound_range_constraints_impl!(
    container::PSI.OptimizationContainer,
    ::Type{T},
    array,
    ::Type{P},
    devices::IS.FlattenIteratorWrapper{V},
    model::DeviceModel{V, W},
) where {
    T <: PSI.ConstraintType,
    P <: PSI.ParameterType,
    V <: PSY.Component,
    W <: FixDispatchablePowerLoad,
}
    time_steps = PSI.get_time_steps(container)
    names = [PSY.get_name(d) for d in devices]

    constraint = PSI.add_constraints_container!(container, T(), V, names, time_steps, meta="ub")

    parameter = PSI.get_parameter_array(container, P(), V)
    multiplier = PSI.get_parameter_multiplier_array(container, P(), V)
    jump_model = PSI.get_jump_model(container)
    fixed_bound_range_with_parameter!(
        jump_model,
        constraint,
        array,
        multiplier,
        parameter,
        devices,
    )
    return
end

function add_parameterized_fixed_bound_range_constraints(
    container::PSI.OptimizationContainer,
    ::Type{T},
    ::Type{U},
    ::Type{P},
    devices::IS.FlattenIteratorWrapper{V},
    model::DeviceModel{V, W},
    ::Type{X},
) where {
    T <: PSI.ConstraintType,
    U <: PSI.VariableType,
    P <: PSI.ParameterType,
    V <: PSY.Component,
    W <: FixDispatchablePowerLoad,
    X <: PM.AbstractPowerModel,
}
    array = PSI.get_variable(container, U(), V)
    _add_parameterized_fixed_bound_range_constraints_impl!(
        container,
        T,
        array,
        P,
        devices,
        model,
    )
    return
end

function PSI.add_constraints!(
    container::PSI.OptimizationContainer,
    T::Type{PSI.ActivePowerVariableLimitsConstraint},
    U::Type{<:PSI.VariableType},
    devices::IS.FlattenIteratorWrapper{V},
    model::DeviceModel{V, W},
    X::Type{<:PM.AbstractPowerModel},
) where {V <: PSY.ControllableLoad, W <: FixDispatchablePowerLoad}
    add_parameterized_fixed_bound_range_constraints(
        container,
        ActivePowerVariableTimeSeriesLimitsConstraint,
        U,
        ActivePowerTimeSeriesParameter,
        devices,
        model,
        X,
    )
    return
end


function fixed_bound_range_with_parameter_60perc!(
    jump_model::JuMP.Model,
    constraint_container::PSI.JuMPConstraintArray,
    lhs_array,
    param_multiplier::PSI.JuMPFloatArray,
    param_array::Union{PSI.JuMPParamArray, PSI.JuMPFloatArray},
    devices::IS.FlattenIteratorWrapper{V},
) where {V <: PSY.Component}
    time_steps = axes(constraint_container)[2]
    for device in devices, t in time_steps
        name = PSY.get_name(device)
        constraint_container[name, t] = JuMP.@constraint(
            jump_model,
            lhs_array[name, t] == 0.6 * param_multiplier[name, t] * param_array[name, t]
        )
    end
    return
end

function _add_parameterized_fixed_bound_range_constraints_impl_60perc!(
    container::PSI.OptimizationContainer,
    ::Type{T},
    array,
    ::Type{P},
    devices::IS.FlattenIteratorWrapper{V},
    model::DeviceModel{V, W},
) where {
    T <: PSI.ConstraintType,
    P <: PSI.ParameterType,
    V <: PSY.Component,
    W <: FixDispatchablePowerLoad_60perc,
}
    time_steps = PSI.get_time_steps(container)
    names = [PSY.get_name(d) for d in devices]

    constraint = PSI.add_constraints_container!(container, T(), V, names, time_steps, meta="ub")

    parameter = PSI.get_parameter_array(container, P(), V)
    multiplier = PSI.get_parameter_multiplier_array(container, P(), V)
    jump_model = PSI.get_jump_model(container)
    fixed_bound_range_with_parameter_60perc!(
        jump_model,
        constraint,
        array,
        multiplier,
        parameter,
        devices,
    )
    return
end

function add_parameterized_fixed_bound_range_constraints_60perc(
    container::PSI.OptimizationContainer,
    ::Type{T},
    ::Type{U},
    ::Type{P},
    devices::IS.FlattenIteratorWrapper{V},
    model::DeviceModel{V, W},
    ::Type{X},
) where {
    T <: PSI.ConstraintType,
    U <: PSI.VariableType,
    P <: PSI.ParameterType,
    V <: PSY.Component,
    W <: FixDispatchablePowerLoad_60perc,
    X <: PM.AbstractPowerModel,
}
    array = PSI.get_variable(container, U(), V)
    _add_parameterized_fixed_bound_range_constraints_impl_60perc!(
        container,
        T,
        array,
        P,
        devices,
        model,
    )
    return
end

function PSI.add_constraints!(
    container::PSI.OptimizationContainer,
    T::Type{PSI.ActivePowerVariableLimitsConstraint},
    U::Type{<:PSI.VariableType},
    devices::IS.FlattenIteratorWrapper{V},
    model::DeviceModel{V, W},
    X::Type{<:PM.AbstractPowerModel},
) where {V <: PSY.ControllableLoad, W <: FixDispatchablePowerLoad_60perc}
    add_parameterized_fixed_bound_range_constraints_60perc(
        container,
        ActivePowerVariableTimeSeriesLimitsConstraint,
        U,
        ActivePowerTimeSeriesParameter,
        devices,
        model,
        X,
    )
    return
end

PSI.sos_status(::FlexibleLoad, ::PSI.AbstractControllablePowerLoadFormulation)=PSI.SOSStatusVariable.NO_VARIABLE
PSI.uses_compact_power(::FlexibleLoad, ::PSI.AbstractControllablePowerLoadFormulation)=false
PSI.sos_status(::FixedFlexibleLoad, ::PSI.AbstractControllablePowerLoadFormulation)=PSI.SOSStatusVariable.NO_VARIABLE
PSI.uses_compact_power(::FixedFlexibleLoad, ::PSI.AbstractControllablePowerLoadFormulation)=false
PSI.sos_status(::CurtailableLoad, ::PSI.AbstractControllablePowerLoadFormulation)=PSI.SOSStatusVariable.NO_VARIABLE
PSI.uses_compact_power(::CurtailableLoad, ::PSI.AbstractControllablePowerLoadFormulation)=false
PSI.sos_status(::ChlorAlkaliPlantWithBidLoad, ::PSI.AbstractControllablePowerLoadFormulation)=PSI.SOSStatusVariable.NO_VARIABLE
PSI.uses_compact_power(::ChlorAlkaliPlantWithBidLoad, ::PSI.AbstractControllablePowerLoadFormulation)=false

function PSI.objective_function!(
    container::PSI.OptimizationContainer,
    devices::IS.FlattenIteratorWrapper{T},
    ::DeviceModel{T, U},
    ::Type{<:PM.AbstractPowerModel},
) where {T <: PSY.ControllableLoad, U <: PSI.AbstractControllablePowerLoadFormulation}
    PSI.add_variable_cost!(container, ActivePowerVariable(), devices, U())
    return
end


function PSI._add_pwl_term!(
    container::PSI.OptimizationContainer,
    component::T,
    cost_data::Matrix{PSY.VariableCost{Vector{Tuple{Float64, Float64}}}},
    ::U,
    ::V,
) where {T <: Union{FlexibleLoad, FixedFlexibleLoad, CurtailableLoad, ChlorAlkaliPlantWithBidLoad}, U <: PSI.VariableType, V <: PSI.AbstractControllablePowerLoadFormulation}
    multiplier = PSI.objective_function_multiplier(U(), V())
    resolution = PSI.get_resolution(container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    base_power = PSI.get_base_power(container)
    # Re-scale breakpoints by Basepower
    name = PSY.get_name(component)
    time_steps = PSI.get_time_steps(container)
    pwl_cost_expressions = Vector{JuMP.AffExpr}(undef, time_steps[end])
    sos_val = PSI._get_sos_value(container, V, component)
    for t in time_steps
        data = PSY.get_cost(cost_data[t])
        slopes = PSY.get_slopes(data)
        # First element of the return is the average cost at P_min.
        # Shouldn't be passed for convexity check
        is_convex = PSI._slope_convexity_check(slopes[2:end])
        break_points = map(x -> last(x), data) ./ base_power
        PSI._add_pwl_variables!(container, T, name, t, data)
        PSI._add_pwl_constraint!(container, component, U(), break_points, sos_val, t)
        if !is_convex
            PSI._add_pwl_sos_constraint!(container, component, U(), break_points, sos_val, t)
        end
        pwl_cost = PSI._get_pwl_cost_expression(container, component, t, data, multiplier * dt)
        pwl_cost_expressions[t] = pwl_cost
    end
    return pwl_cost_expressions
end


function PSI.add_parameters!(
    container::PSI.OptimizationContainer,
    ::T,
    devices::U,
    model::PSI.DeviceModel{D, W},
) where {
    T <: PSI.TimeSeriesParameter,
    U <: Union{Vector{D}, IS.FlattenIteratorWrapper{D}},
    W <: PSI.AbstractControllablePowerLoadFormulation,
} where {D <: Union{FlexibleLoad, FixedFlexibleLoad, CurtailableLoad, ChlorAlkaliPlantWithBidLoad}}
    ts_type = PSI.get_default_time_series_type(container)
    if !(ts_type <: Union{PSY.AbstractDeterministic, PSY.StaticTimeSeries})
        error("add_parameters! for TimeSeriesParameter is not compatible with $ts_type")
    end
    time_steps = PSI.get_time_steps(container)
    names = [PSY.get_name(d) for d in devices]
    ts_name = PSI.get_time_series_names(model)[T]
    time_series_mult_id = PSI.create_time_series_multiplier_index(model, T)
    @debug "adding" T ts_name ts_type time_series_mult_id _group =
        :OptimizationContainer
    parameter_container =
        PSI.add_param_container!(container, T(), D, ts_type, ts_name, names, time_steps)
    PSI.set_time_series_multiplier_id!(PSI.get_attributes(parameter_container), time_series_mult_id)
    jump_model = PSI.get_jump_model(container)
    for d in devices
        name = PSY.get_name(d)
        ts_vector = PSI.get_time_series(container, d, T())
        multiplier = PSI.get_multiplier_value(T(), d, W())
        if get_name(d) == "Clark"
            @show ts_vector
            @show multiplier
        end
        for t in time_steps
            PSI.set_parameter!(
                parameter_container,
                jump_model,
                ts_vector[t],
                multiplier,
                name,
                t,
            )
        end
    end
    return
end


function PSI._update_parameter_values!(
    param_array,
    attributes::PSI.CostFunctionAttributes,
    ::Type{V},
    model::PSI.DecisionModel,
    ::PSI.DatasetContainer{PSI.DataFrameDataset},
) where {V <: Union{FlexibleLoad, FixedFlexibleLoad, CurtailableLoad, ChlorAlkaliPlantWithBidLoad}}
    initial_forecast_time = PSI.get_current_time(model) # Function not well defined for DecisionModels
    time_steps = PSI.get_time_steps(PSI.get_optimization_container(model))
    horizon = time_steps[end]
    container = PSI.get_optimization_container(model)
    if PSI.is_synchronized(container)
        obj_func = PSI.get_objective_function(container)
        PSI.set_synchronized_status(obj_func, false)
        PSI.reset_variant_terms(obj_func)
    end
    components = PSI.get_available_components(V, PSI.get_system(model))

    for component in components
        if PSI._has_variable_cost_parameter(component)
            name = PSY.get_name(component)
            ts_vector = PSY.get_variable_cost(
                component,
                PSY.get_operation_cost(component);
                start_time=initial_forecast_time,
                len=horizon,
            )
            variable_cost_forecast_values = TimeSeries.values(ts_vector)
            for (t, value) in enumerate(variable_cost_forecast_values)
                if attributes.uses_compact_power
                    value, _ = PSI._convert_variable_cost(value)
                end
                PSI._set_param_value!(param_array, PSY.get_cost(value), name, t)
                PSI.update_variable_cost!(container, param_array, attributes, component, t)
            end
        end
    end

    jump_model = PSI.get_jump_model(container)

    if V <: FlexibleLoad
        for pwl_constraint in PSI.get_constraints(container)[PSI.ConstraintKey{PSI.PieceWiseLinearCostConstraint, FlexibleLoad}("")]
            JuMP.delete(jump_model, pwl_constraint)
        end
    elseif V <: FixedFlexibleLoad
        for pwl_constraint in PSI.get_constraints(container)[PSI.ConstraintKey{PSI.PieceWiseLinearCostConstraint, FixedFlexibleLoad}("")]
            JuMP.delete(jump_model, pwl_constraint)
        end
    elseif V <: CurtailableLoad
        for pwl_constraint in PSI.get_constraints(container)[PSI.ConstraintKey{PSI.PieceWiseLinearCostConstraint, CurtailableLoad}("")]
            JuMP.delete(jump_model, pwl_constraint)
        end
    elseif V <: ChlorAlkaliPlantWithBidLoad
        for pwl_constraint in PSI.get_constraints(container)[PSI.ConstraintKey{PSI.PieceWiseLinearCostConstraint, ChlorAlkaliPlantWithBidLoad}("")]
            JuMP.delete(jump_model, pwl_constraint)
        end
    end

    pwl_vars = PSI.get_variable(container, PSI.PieceWiseLinearCostVariable(), V)
    variables = PSI.get_variable(container, PSI.ActivePowerVariable(), V)
    const_container = PSI.get_constraint(
            container,
            PSI.PieceWiseLinearCostConstraint(),
            V,
    )
    base_power = PSI.get_base_power(container)

    for component in components
        name = PSY.get_name(component)
        ts_vector = PSY.get_variable_cost(
                component,
                PSY.get_operation_cost(component);
                start_time=initial_forecast_time,
                len=horizon,
            )
        variable_cost_forecast_values = TimeSeries.values(ts_vector)

        for (t, value) in enumerate(variable_cost_forecast_values)
            len_cost_data = length(value)
            break_point = Float64[]
            for i in 1:len_cost_data
                push!(break_point, value[i][2] / base_power)
            end

            const_container[name, t] = JuMP.@constraint(
                jump_model,
                variables[name, t] ==
                sum(pwl_vars[name, ix, t] * break_point[ix] for ix in 1:len_cost_data)
            )
        end
    end

    return
end


function PSI._update_pwl_cost_expression(
    container::PSI.OptimizationContainer,
    ::Type{T},
    component_name::String,
    time_period::Int,
    cost_data::Vector{NTuple{2, Float64}},
) where {T <: Union{FlexibleLoad, FixedFlexibleLoad, CurtailableLoad, ChlorAlkaliPlantWithBidLoad}}
    pwl_var_container = PSI.get_variable(container, PSI.PieceWiseLinearCostVariable(), T)
    resolution = PSI.get_resolution(container)
    dt = Dates.value(Dates.Second(resolution)) / SECONDS_IN_HOUR
    gen_cost = JuMP.AffExpr(0.0)
    for i in 1:length(cost_data)
        JuMP.add_to_expression!(
            gen_cost,
            cost_data[i][1] * dt * pwl_var_container[(component_name, i, time_period)],
        )
    end
    return gen_cost
end


function PSI.add_constraints!(
    container::PSI.OptimizationContainer,
    T::Type{PSI.ActivePowerVariableLimitsConstraint},
    U::Type{<:PSI.VariableType},
    devices::IS.FlattenIteratorWrapper{V},
    model::PSI.DeviceModel{V, W},
    X::Type{<:PM.AbstractPowerModel},
) where {V <: PSY.ControllableLoad, W <: CurtailablePowerLoad}
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

function PSI.objective_function!(
    container::PSI.OptimizationContainer,
    devices::IS.FlattenIteratorWrapper{T},
    ::PSI.DeviceModel{T, U},
    ::Type{<:PM.AbstractPowerModel},
) where {T <: PSY.ControllableLoad, U <: CurtailablePowerLoad}
    PSI.add_variable_cost!(container, PSI.ActivePowerVariable(), devices, U())
    return
end


PSI.get_variable_multiplier(_, ::Type{CurtailableLoad}, ::CurtailablePowerLoad) = 1.0

function PSI.get_default_time_series_names(
    ::Type{<:ChlorAlkaliPlantWithBidLoad},
    ::Type{<:PSI.AbstractControllablePowerLoadFormulation},
)
    return Dict{Type{<:PSI.TimeSeriesParameter}, String}(
        PSI.ActivePowerTimeSeriesParameter => "max_active_power",
        PSI.EnergyBudgetTimeSeriesParameter => "load_budget",
    )
end


function PSI.construct_device!(
    container::PSI.OptimizationContainer,
    sys::PSY.System,
    ::PSI.ArgumentConstructStage,
    model::PSI.DeviceModel{L, D},
    ::Type{S},
) where {
    L <: ChlorAlkaliPlantWithBidLoad,
    D <: ChlorAlkaliPlantDispatchWithBid,
    S <: PM.AbstractActivePowerModel,
}
    devices = PSI.get_available_components(L, sys)

    PSI.add_variables!(container, ActivePowerVariable, devices, D())
    PSI.add_parameters!(container, PSI.EnergyBudgetTimeSeriesParameter, devices, model)

    # Add Variables to expressions
    PSI.add_to_expression!(
        container,
        ActivePowerBalance,
        ActivePowerVariable,
        devices,
        model,
        S,
    )

    PSI.add_parameters!(container, PSI.ActivePowerTimeSeriesParameter, devices, model)

    PSI.add_expressions!(container, PSI.ProductionCostExpression, devices, model)
    return
end

function PSI.construct_device!(
    container::PSI.OptimizationContainer,
    sys::PSY.System,
    ::PSI.ModelConstructStage,
    model::PSI.DeviceModel{L, D},
    ::Type{S},
) where {
    L <: ChlorAlkaliPlantWithBidLoad,
    D <: ChlorAlkaliPlantDispatchWithBid,
    S <: PM.AbstractActivePowerModel,
}
    devices = PSI.get_available_components(L, sys)

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

    PSI.objective_function!(container, devices, model, S)

    PSI.add_constraint_dual!(container, sys, model)
    return
end


function PSI.get_min_max_limits(
    x::ChlorAlkaliPlantWithBidLoad,
    ::Type{<:PSI.ActivePowerVariableLimitsConstraint},
    ::Type{ChlorAlkaliPlantDispatchWithBid},
)
    return PSY.get_active_power_limits(x)
end