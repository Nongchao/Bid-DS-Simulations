# ### Modeling Packages
using PowerSystems
const PSY = PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using PowerSystemCaseBuilder
using InfrastructureSystems
const IS = InfrastructureSystems
using PowerModels
const PM = PowerModels

# ### Data management packages
using Dates
using DataFrames
using TimeSeries
using CSV
using DelimitedFiles
using Statistics

# ### Optimization packages
using Xpress
using JuMP
using Gurobi
# const GRB_ENV = Gurobi.Env()

# ### Plotting packages
# using PlotlyJS
# using Plots

# ### Other packages
using BenchmarkTools
using Logging
using Pkg

# ### Knobs to tweak
carbon_pricing = "false"   #ARGS[1]
load_model = "fixed_load"   # "fixed_load", "google_CAC", "iso_load_shift"
peak_load = "capped"    # "capped", "uncapped"
number_of_data_center = "All"  # "Barton", "All"
days = 365
horizon = 24
interval = Dates.Hour(24)
start_date = Dates.DateTime("2020-01-01T00:00:00")
fixed_load_perc = 0.5
MIP_solver = Xpress
Google_emission_cost = 0.02268
Google_capacity_cost = 50.0

proj = Pkg.project().path
proj_path = dirname(proj)

result_file_path = joinpath(proj_path, "results", "google_model_"*"$load_model"*"_carbon_price_"*"$carbon_pricing"*"_steps_"*"$days"*"_updated_"*"$peak_load"*"_data_center_"*"$number_of_data_center")
mkpath(result_file_path)
scratch_output = "/lustre/eaglefs/scratch/nguo/project/google_output"
mkpath(scratch_output)


# ### functions to include
include(joinpath(proj_path, "src", "bid_model_function.jl"))
include(joinpath(proj_path, "src", "shiftable_load_function.jl"))
include(joinpath(proj_path, "src", "Google_CAC_function.jl"))
include(joinpath(proj_path, "src", "parse_system.jl"))
include(joinpath(proj_path, "src", "metrics.jl"))

repo = joinpath(proj_path, "..", "RTS-GMLC")
source_gen_data = CSV.read(joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), DataFrame)
if !isdir(joinpath(repo, "gen_original"))
    mkdir(joinpath(repo, "gen_original"))
    cp(joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), joinpath(repo, "gen_original", "gen.csv"), force=true)
end
original_source_gen_data = joinpath(repo, "gen_original", "gen.csv")

emission_type = ["SO2", "NOX", "Part", "CO2", "CH4", "N2O", "CO", "VOCs"]

##### unit: dollars/lb #####
emission_cost = DataFrame()
for e in emission_type
    emission_cost[!, Symbol(e)] = zeros(1)
end
emission_cost[1, :SO2] = 36.5
emission_cost[1, :NOX] = 3.2
emission_cost[1, :Part] = 68.5
# social cost of CO2 is based on $50 per metric ton
# 1 metric ton = 2204.62 lbs
# https://www.whitehouse.gov/wp-content/uploads/2021/02/TechnicalSupportDocument_SocialCostofCarbonMethaneNitrousOxide.pdf?source=email
emission_cost[1, :CO2] = 0.02268

if carbon_pricing == "true"
    add_emission_cost(emission_cost, source_gen_data)
    sys = create_rts_gmlc_system(repo)
    cp(original_source_gen_data, joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), force=true)
else
    sys = create_rts_gmlc_system(repo)
end

for i in get_components(ThermalStandard, sys, x -> occursin("SYNC_COND", x.name))
    remove_component!(sys, i)
end

PTDF_matrix = PTDF(sys)
if number_of_data_center == "Barton"
    data_center = [get_component(PowerLoad, sys, "Barton")]
elseif number_of_data_center == "All"
    data_center = collect(get_components(PowerLoad, sys))
end
non_data_center = collect(get_components(PowerLoad, sys))
for i in data_center
    filter!(x -> x != i, non_data_center)
end

for l in collect(get_components(Line, sys))
    if get_name(l) != "B24"
        set_rate!(l, get_rate(l) * 1)
    else
        set_rate!(l, get_rate(l) * 1)
    end
end

############################################ Fixed Load ##########################################

set_units_base_system!(sys, "NATURAL_UNITS")
remove_time_series!(sys, DeterministicSingleTimeSeries)

set_units_base_system!(sys, "SYSTEM_BASE")

transform_single_time_series!(sys, horizon, interval)

template = ProblemTemplate(
    NetworkModel(
        PSI.StandardPTDFModel,
        PTDF = PTDF(sys),
        duals = [CopperPlateBalanceConstraint],
        use_slacks = true,
    ),
)
set_device_model!(template, DeviceModel(Line, StaticBranch))
set_device_model!(template, DeviceModel(Transformer2W, StaticBranch))
set_device_model!(template, DeviceModel(TapTransformer, StaticBranch))
# set_device_model!(template, DeviceModel(ThermalStandard, ThermalStandardUnitCommitment))
set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicUnitCommitment))
set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
set_device_model!(template, DeviceModel(HydroDispatch, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(HydroEnergyReservoir, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(RenewableFix, FixedOutput))
set_device_model!(template, DeviceModel(PowerLoad, StaticPowerLoad))

# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R1", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R2", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R3", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Reg_Up", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Flex_Up", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, "Reg_Down", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, "Flex_Down", use_slacks = true, duals = [RequirementConstraint]))

for (k, v) in template.branches
    v.duals = [NetworkFlowConstraint]
end

solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.00000001, "OUTPUTLOG" => 1, "PRESOLVE" => 0)
# solver = Gurobi.Optimizer

models = SimulationModels(
    decision_models = [
        DecisionModel(
            template,
            sys,
            optimizer = solver,
            name = "UCED",
            optimizer_solve_log_print = true,
            initialize_model = false,
        ),
    ],
)

sequence = SimulationSequence(
    models = models,
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(
    name = "Bid-DS-RTS-Google",
    steps = days,
    models = models,
    sequence = sequence,
    simulation_folder = scratch_output,
    initial_time = start_date
)

build!(
    sim,
    console_level = Logging.Error,
    file_level = Logging.Debug,
    recorders = [:simulation],
    serialize = false,
)

status = execute!(sim, enable_progress_bar = false, disable_timer_outputs = true)

# f = open("fixed_load_step_$days.txt","w"); print(f, sim.models.decision_models[1].internal.container.JuMPmodel); close(f)

res = get_decision_problem_results(SimulationResults(sim), "UCED")
average_emission_rate_by_region, time_step = results_output(res, sim, PTDF_matrix, emission_type, source_gen_data, result_file_path, start_date, "constant_load")

############################################ Google CAC ##########################################
source_gen_data = CSV.read(joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), DataFrame)

load_model = "google_CAC"   # "fixed_load", "google_CAC", "iso_load_shift"
result_file_path = joinpath("results", "google_model_"*"$load_model"*"_carbon_price_"*"$carbon_pricing"*"_steps_"*"$days"*"_updated_"*"$peak_load"*"_data_center_"*"$number_of_data_center")
mkpath(result_file_path)

if carbon_pricing == "true"
    add_emission_cost(emission_cost, source_gen_data)
    sys = create_rts_gmlc_system(repo)
    cp(original_source_gen_data, joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), force=true)
else
    sys = create_rts_gmlc_system(repo)
end

for i in get_components(ThermalStandard, sys, x -> occursin("SYNC_COND", x.name))
    remove_component!(sys, i)
end

PTDF_matrix = PTDF(sys)
if number_of_data_center == "Barton"
    data_center = [get_component(PowerLoad, sys, "Barton")]
elseif number_of_data_center == "All"
    data_center = collect(get_components(PowerLoad, sys))
end
non_data_center = collect(get_components(PowerLoad, sys))
for i in data_center
    filter!(x -> x != i, non_data_center)
end

for l in collect(get_components(Line, sys))
    if get_name(l) != "B24"
        set_rate!(l, get_rate(l) * 1)
    else
        set_rate!(l, get_rate(l) * 1)
    end
end

set_units_base_system!(sys, "NATURAL_UNITS")
remove_time_series!(sys, DeterministicSingleTimeSeries)

y_all = DataFrame()
for l in data_center
    y_all[!, Symbol(get_name(l))] = Vector{Float64}(undef, sim.steps)
end

for l in data_center
    new_load = zeros(0)
    y_vec = zeros(0)
    region = get_number(get_bus(l))
    for s in 1:sim.steps
        if region < 200
            emission_rate = average_emission_rate_by_region[(s-1)*24+1:s*24, "Region_1"]
        elseif region >= 300
            emission_rate = average_emission_rate_by_region[(s-1)*24+1:s*24, "Region_3"]
        else
            emission_rate = average_emission_rate_by_region[(s-1)*24+1:s*24, "Region_2"]
        end
        
        data_center_load = values(get_time_series_array(SingleTimeSeries, l, "max_active_power"))[(s-1)*24+1:s*24]
        flex_load = (1 - fixed_load_perc) * data_center_load
        fixed_load = data_center_load .- flex_load
        if peak_load == "capped"
            p, y = Google_CAC(emission_rate, Google_emission_cost, Google_capacity_cost, flex_load, fixed_load, MIP_solver)
            append!(y_vec, y)
        else
            p = Google_CAC_No_Peak(emission_rate, Google_emission_cost, get_max_active_power(l), flex_load, fixed_load, MIP_solver)
        end
        Google_CAC_load = p .+ fixed_load
        # append!(y_vec, y)
        append!(new_load, Google_CAC_load)
    end
    if peak_load == "capped"
        y_all[!, Symbol(get_name(l))] = y_vec
    end
    timestamp = collect(DateTime("2020-01-01T00:00:00"):Hour(1):(DateTime("2020-01-01T00:00:00") + Hour(8783)));
    load_ts = values(get_time_series_array(SingleTimeSeries, l, "max_active_power"))
    load_ts[1:24*sim.steps] = new_load
    remove_time_series!(sys, SingleTimeSeries, l, "max_active_power")

    load_data = TimeArray(timestamp, load_ts ./ get_max_active_power(l))
    _load_data_time_series = SingleTimeSeries("max_active_power", load_data, scaling_factor_multiplier=get_max_active_power)
    add_time_series!(sys, l, _load_data_time_series)
end

set_units_base_system!(sys, "SYSTEM_BASE")

transform_single_time_series!(sys, horizon, interval)

template = ProblemTemplate(
    NetworkModel(
        PSI.StandardPTDFModel,
        PTDF = PTDF(sys),
        duals = [CopperPlateBalanceConstraint],
        use_slacks = true,
    ),
)
set_device_model!(template, DeviceModel(Line, StaticBranch))
set_device_model!(template, DeviceModel(Transformer2W, StaticBranch))
set_device_model!(template, DeviceModel(TapTransformer, StaticBranch))
# set_device_model!(template, DeviceModel(ThermalStandard, ThermalStandardUnitCommitment))
set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicUnitCommitment))
set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
set_device_model!(template, DeviceModel(HydroDispatch, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(HydroEnergyReservoir, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(RenewableFix, FixedOutput))
set_device_model!(template, DeviceModel(PowerLoad, StaticPowerLoad))

# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R1", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R2", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R3", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Reg_Up", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Flex_Up", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, "Reg_Down", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, "Flex_Down", use_slacks = true, duals = [RequirementConstraint]))

for (k, v) in template.branches
    v.duals = [NetworkFlowConstraint]
end

solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.00000001, "OUTPUTLOG" => 1, "PRESOLVE" => 0)
# solver = Gurobi.Optimizer

models = SimulationModels(
    decision_models = [
        DecisionModel(
            template,
            sys,
            optimizer = solver,
            name = "UCED",
            optimizer_solve_log_print = true,
            initialize_model = false,
        ),
    ],
)

sequence = SimulationSequence(
    models = models,
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(
    name = "Bid-DS-RTS-Google",
    steps = days,
    models = models,
    sequence = sequence,
    simulation_folder = scratch_output,
    initial_time = start_date
)

build!(
    sim,
    console_level = Logging.Error,
    file_level = Logging.Debug,
    recorders = [:simulation],
    serialize = false,
)

status = execute!(sim, enable_progress_bar = false, disable_timer_outputs = true)

# f = open("google_load_error_233_step_$days.txt","w"); print(f, sim.models.decision_models[1].internal.container.JuMPmodel); close(f)

res = get_decision_problem_results(SimulationResults(sim), "UCED")
results_output(res, sim, PTDF_matrix, emission_type, source_gen_data, result_file_path, start_date, "constant_load")

############################################ ISO Dispatch ##########################################
source_gen_data = CSV.read(joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), DataFrame)

load_model = "iso_load_shift"   # "fixed_load", "google_CAC", "iso_load_shift"
result_file_path = joinpath("results", "google_model_"*"$load_model"*"_carbon_price_"*"$carbon_pricing"*"_steps_"*"$days"*"_updated_"*"$peak_load"*"_data_center_"*"$number_of_data_center")
mkpath(result_file_path)

if carbon_pricing == "true"
    add_emission_cost(emission_cost, source_gen_data)
    sys = create_rts_gmlc_system(repo)
    cp(original_source_gen_data, joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), force=true)
else
    sys = create_rts_gmlc_system(repo)
end

for i in get_components(ThermalStandard, sys, x -> occursin("SYNC_COND", x.name))
    remove_component!(sys, i)
end

PTDF_matrix = PTDF(sys)
if number_of_data_center == "Barton"
    data_center = [get_component(PowerLoad, sys, "Barton")]
elseif number_of_data_center == "All"
    data_center = collect(get_components(PowerLoad, sys))
end
non_data_center = collect(get_components(PowerLoad, sys))
for i in data_center
    filter!(x -> x != i, non_data_center)
end

for l in collect(get_components(Line, sys))
    if get_name(l) != "B24"
        set_rate!(l, get_rate(l) * 1)
    else
        set_rate!(l, get_rate(l) * 1)
    end
end

set_units_base_system!(sys, "NATURAL_UNITS")
remove_time_series!(sys, DeterministicSingleTimeSeries)

# use peak load as upper bound
if peak_load == "uncapped"
    y_all = DataFrame()
    for l in data_center
        y_all[!, Symbol(get_name(l))] = Float64[]
    end
    y_all_dict = Dict()
    for l in data_center
        push!(y_all_dict, Symbol(get_name(l)) => get_max_active_power(l))       
    end
    for s in 1:sim.steps
        push!(y_all, y_all_dict)
    end
end

# use peak load as upper bound
# for l in data_center
#     y_all[!, Symbol(get_name(l))] = ones(sim.steps) * get_max_active_power(l)
# end

# replace PowerLoad with ShiftableLoad
for l in data_center
    # construct bid time series
    timestamp = collect(DateTime("2020-01-01T00:00:00"):Hour(1):(DateTime("2020-01-01T00:00:00") + Hour(8783)));
    load_ts = values(get_time_series_array(SingleTimeSeries, l, "max_active_power")) * fixed_load_perc
    load_ts_flex = values(get_time_series_array(SingleTimeSeries, l, "max_active_power")) * (1 - fixed_load_perc)
    remove_time_series!(sys, SingleTimeSeries, l, "max_active_power")

    load_data = TimeArray(timestamp, load_ts ./ get_max_active_power(l))
    _load_data_time_series = SingleTimeSeries("max_active_power", load_data, scaling_factor_multiplier=get_max_active_power)
    add_time_series!(sys, l, _load_data_time_series)

    # replace PowerLoad with ShiftableLoad
    base_power = get_base_power(l)
    load_budget_data = TimeArray(timestamp, load_ts_flex ./ get_max_active_power(l))
    _load_budget_time_series = SingleTimeSeries("load_budget", load_budget_data, scaling_factor_multiplier=get_max_active_power)

    max_hourly_load_budget_data_value = copy(load_ts)
    for s in 1:sim.steps
        max_hourly_load_budget_data_value[(s-1)*24+1:s*24] = y_all[s, Symbol(get_name(l))] .- load_ts[(s-1)*24+1:s*24]
    end
    max_hourly_load_budget_data = TimeArray(timestamp, max_hourly_load_budget_data_value ./ get_max_active_power(l))
    _max_hourly_load_budget_time_series = SingleTimeSeries("max_active_power", max_hourly_load_budget_data, scaling_factor_multiplier=get_max_active_power)

    flex_load = ShiftableLoad(
        name = get_name(l),
        available = true,
        bus = get_bus(l),
        model = get_model(l),   # Is this right?
        active_power = get_active_power(l) / base_power,
        reactive_power = get_reactive_power(l) / base_power,
        max_active_power = get_active_power(l) / base_power,
        max_reactive_power = get_reactive_power(l) / base_power,
        active_power_limits = (min = 0.0, max = get_active_power(l) / base_power),
        base_power = base_power,
        services = get_services(l),
        dynamic_injector = get_dynamic_injector(l),
        ext = get_ext(l),
    )
    add_component!(sys, flex_load)
    add_time_series!(sys, flex_load, _load_budget_time_series)
    add_time_series!(sys, flex_load, _max_hourly_load_budget_time_series)
    # we do not need to set available to false because the original fixed part is still there
end

set_units_base_system!(sys, "SYSTEM_BASE")

transform_single_time_series!(sys, horizon, interval)

template = ProblemTemplate(
    NetworkModel(
        PSI.StandardPTDFModel,
        PTDF = PTDF(sys),
        duals = [CopperPlateBalanceConstraint],
        use_slacks = true,
    ),
)
set_device_model!(template, DeviceModel(Line, StaticBranch))
set_device_model!(template, DeviceModel(Transformer2W, StaticBranch))
set_device_model!(template, DeviceModel(TapTransformer, StaticBranch))
# set_device_model!(template, DeviceModel(ThermalStandard, ThermalStandardUnitCommitment))
set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicUnitCommitment))
set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
set_device_model!(template, DeviceModel(HydroDispatch, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(HydroEnergyReservoir, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(RenewableFix, FixedOutput))
set_device_model!(template, DeviceModel(PowerLoad, StaticPowerLoad))

set_device_model!(template, DeviceModel(ShiftableLoad, LoadShifting))

# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R1", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R2", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Spin_Up_R3", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Reg_Up", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, "Flex_Up", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, "Reg_Down", use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, "Flex_Down", use_slacks = true, duals = [RequirementConstraint]))

for (k, v) in template.branches
    v.duals = [NetworkFlowConstraint]
end

solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.00000001, "OUTPUTLOG" => 1, "PRESOLVE" => 0)
# solver = Gurobi.Optimizer

models = SimulationModels(
    decision_models = [
        DecisionModel(
            template,
            sys,
            optimizer = solver,
            name = "UCED",
            optimizer_solve_log_print = true,
            initialize_model = false,
        ),
    ],
)

sequence = SimulationSequence(
    models = models,
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(
    name = "Bid-DS-RTS-Google",
    steps = days,
    models = models,
    sequence = sequence,
    simulation_folder = scratch_output,
    initial_time = start_date
)

build!(
    sim,
    console_level = Logging.Error,
    file_level = Logging.Debug,
    recorders = [:simulation],
    serialize = false,
)

status = execute!(sim, enable_progress_bar = false, disable_timer_outputs = true)

f = open("iso_shift_step_$days.txt","w"); print(f, sim.models.decision_models[1].internal.container.JuMPmodel); close(f)

res = get_decision_problem_results(SimulationResults(sim), "UCED")
results_output(res, sim, PTDF_matrix, emission_type, source_gen_data, result_file_path, start_date, "shiftable_load")