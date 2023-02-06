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
const GRB_ENV = Gurobi.Env()
const Gurobi_optimizer = optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV))

# ### Plotting packages
# using PlotlyJS
# using Plots

# ### Other packages
using BenchmarkTools
using Logging
using Pkg

# ### Knobs to tweak
carbon_pricing = "false"  #ARGS[1]
congestion = "original_congestion" # "no_congestion", "with_congestion", "original_congestion"
number_of_data_center = "Barton"  # "Barton", "All"
days = 2
horizon = 24
interval = Dates.Hour(24)
start_date = Dates.DateTime("2020-01-01T00:00:00")
thermal_formulation = "basic" # "basic", "standard"
const chlor_alkali_bids = 330.0
const chlor_alkali_daily_consumption = 55.0 / 0.422
const chlor_alkali_max_load = 9.0
const chlor_alkali_min_load = 3.77

proj = Pkg.project().path
proj_path = dirname(proj)

result_file_path = joinpath(proj_path, "results", "chlor_alkalai_bid_model"*"_$congestion"*"_carbon_price_"*"$carbon_pricing"*"_steps_"*"$days"*"_thermal_formulation_"*"$thermal_formulation"*"_data_center_"*"$number_of_data_center")
mkpath(result_file_path)
scratch_output = "/lustre/eaglefs/scratch/nguo/project/google_output"
mkpath(scratch_output)

# ### functions to include
include(joinpath(proj_path, "src", "bid_model_function.jl"))
include(joinpath(proj_path, "src", "shiftable_load_function.jl"))
include(joinpath(proj_path, "src", "parse_system.jl"))
include(joinpath(proj_path, "src", "metrics_chlor_alkali.jl"))

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

if congestion == "no_congestion"
    for l in collect(get_components(Line, sys))
        set_rate!(l, get_rate(l) * 100)
    end
elseif congestion == "with_congestion"
    for l in collect(get_components(Line, sys))
        if get_name(l) != "B24"
            set_rate!(l, get_rate(l) * 100)
        else
            set_rate!(l, get_rate(l) * 0.1)
        end
    end
end

set_units_base_system!(sys, "NATURAL_UNITS")
remove_time_series!(sys, DeterministicSingleTimeSeries)


# for l in data_center # get_components(PowerLoad, sys)
#     base_power = get_base_power(l)
#     # construct bid time series
#     timestamp = collect(DateTime("2020-01-01T00:00:00"):Hour(1):(DateTime("2020-01-01T00:00:00") + Hour(8783)));
#     load_ts = get_time_series_values(SingleTimeSeries, l, "max_active_power")
#     bid_ts = collect([(0.0, 0.0), (-330*load_ts[t], load_ts[t])] for t = 1:length(load_ts))
#     bid_data = TimeArray(timestamp, bid_ts)
#     _bid_time_series = SingleTimeSeries("variable_cost", bid_data)

#     load_ts_flex = values(get_time_series_array(SingleTimeSeries, l, "max_active_power")) * 0.8
#     load_budget_data = TimeArray(timestamp, load_ts_flex ./ get_max_active_power(l))
#     _load_budget_time_series = SingleTimeSeries("load_budget", load_budget_data, scaling_factor_multiplier=get_max_active_power)

#     flex_load = ChlorAlkaliPlantWithBidLoad(
#         name = "Chlor_Alkali_Bid_"*get_name(l),
#         available = true,
#         bus = get_bus(l),
#         model = get_model(l),   # Is this right?
#         active_power = get_active_power(l) / base_power,
#         reactive_power = get_reactive_power(l) / base_power,
#         max_active_power = get_active_power(l) / base_power,
#         max_reactive_power = get_reactive_power(l) / base_power,
#         active_power_limits = (min = 0.0, max = get_active_power(l) / base_power),
#         base_power = base_power,
#         operation_cost = MarketBidCost(0.0, (hot = 0.0, warm = 0.0, cold = 0.0), 0.0, InfrastructureSystems.TimeSeriesKey(_bid_time_series), Vector{Service}[]),   # Is this just going to over-written by the time series? This is the total cost - quantity pair.
#         services = get_services(l),
#         dynamic_injector = get_dynamic_injector(l),
#         ext = get_ext(l),
#     )
#     add_component!(sys, flex_load)
#     copy_time_series!(flex_load, l)
#     add_time_series!(sys, flex_load, _bid_time_series)
#     add_time_series!(sys, flex_load, _load_budget_time_series)
# end



for l in data_center # get_components(PowerLoad, sys)
    base_power = get_base_power(l)
    # construct bid time series
    timestamp = collect(DateTime("2020-01-01T00:00:00"):Hour(1):(DateTime("2020-01-01T00:00:00") + Hour(8783)));
    load_data = TimeArray(timestamp, ones(8784))
    _load_data_time_series = SingleTimeSeries("max_active_power", load_data, scaling_factor_multiplier=get_max_active_power)

    bid_ts = collect([(0.0, 0.0), (- chlor_alkali_bids * chlor_alkali_max_load, chlor_alkali_max_load)] for t = 1:8784)
    bid_data = TimeArray(timestamp, bid_ts)
    _bid_time_series = SingleTimeSeries("variable_cost", bid_data)

    load_ts_flex = ones(8784) * chlor_alkali_daily_consumption / 24
    load_budget_data = TimeArray(timestamp, load_ts_flex ./ chlor_alkali_max_load)
    _load_budget_time_series = SingleTimeSeries("load_budget", load_budget_data, scaling_factor_multiplier=get_max_active_power)

    flex_load = ChlorAlkaliPlantWithBidLoad(
        name = "Chlor_Alkali_Bid_"*get_name(l),
        available = true,
        bus = get_bus(l),
        model = get_model(l),   # Is this right?
        active_power = chlor_alkali_max_load / base_power,
        reactive_power = 0.0,
        max_active_power = chlor_alkali_max_load / base_power,
        max_reactive_power = 0.0,
        active_power_limits = (min = chlor_alkali_min_load / base_power, max = chlor_alkali_max_load / base_power),
        base_power = base_power,
        operation_cost = MarketBidCost(0.0, (hot = 0.0, warm = 0.0, cold = 0.0), 0.0, InfrastructureSystems.TimeSeriesKey(_bid_time_series), Vector{Service}[]),   # Is this just going to over-written by the time series? This is the total cost - quantity pair.
        services = get_services(l),
        dynamic_injector = get_dynamic_injector(l),
        ext = get_ext(l),
    )
    add_component!(sys, flex_load)
    add_time_series!(sys, flex_load, _load_data_time_series)
    add_time_series!(sys, flex_load, _bid_time_series)
    add_time_series!(sys, flex_load, _load_budget_time_series)
end
set_units_base_system!(sys, "SYSTEM_BASE")


transform_single_time_series!(sys, horizon, interval)

template = ProblemTemplate(
    NetworkModel(
        PSI.StandardPTDFModel,
        PTDF = PTDF(sys),
        duals = [CopperPlateBalanceConstraint],
        use_slacks = false,
    ),
)
set_device_model!(template, DeviceModel(Line, StaticBranch))
set_device_model!(template, DeviceModel(TapTransformer, StaticBranch))
if thermal_formulation == "basic"
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicUnitCommitment))
elseif thermal_formulation == "standard"
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalStandardUnitCommitment))
end
set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
set_device_model!(template, DeviceModel(HydroDispatch, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(HydroEnergyReservoir, HydroDispatchRunOfRiver))
set_device_model!(template, DeviceModel(RenewableFix, FixedOutput))
set_device_model!(template, DeviceModel(PowerLoad, StaticPowerLoad))
set_device_model!(template, DeviceModel(ChlorAlkaliPlantWithBidLoad, ChlorAlkaliPlantDispatchWithBid))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve, use_slacks = true, duals = [RequirementConstraint]))
# set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve, use_slacks = true, duals = [RequirementConstraint]))
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

solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.0001, "OUTPUTLOG" => 1, "PRESOLVE" => 0)
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
            export_pwl_vars = true,
        ),
    ],
)

sequence = SimulationSequence(
    models = models,
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(
    name = "Bid-DS-RTS",
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

status = execute!(sim, enable_progress_bar = true, disable_timer_outputs = false)

# f = open("chlor_alkali_formulation_time_series.txt","w"); print(f, sim.models.decision_models[1].internal.container.JuMPmodel); close(f)
# ic = PSI.get_initial_condition(container, PSI.DevicePower(), PSY.ThermalStandard)

res = get_decision_problem_results(SimulationResults(sim), "UCED")
results_output(res, sim, PTDF_matrix, emission_type, source_gen_data, result_file_path, start_date, "bid_in_load", chlor_alkali_bids)