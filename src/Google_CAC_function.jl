function Google_CAC(emission_rate::Vector{Float64}, emission_cost::Float64, capacity_cost::Float64, flex_load::Vector{Float64}, fixed_load::Vector{Float64}, optimizer)
    data_center_model = Model(optimizer.Optimizer)

    @variable(data_center_model, p[1:length(flex_load)] >= 0)
    @variable(data_center_model, y >= 0)

    @objective(data_center_model, Min, emission_cost * sum(emission_rate[h] * p[h] for h in 1:length(flex_load)) + capacity_cost * y)

    @constraint(data_center_model, Eq_total_demand, sum(p[h] for h in 1:length(flex_load)) == sum(flex_load[h] for h in 1:length(flex_load)))
    @constraint(data_center_model, Eq_peak_demand[h in 1:length(flex_load)], p[h] + fixed_load[h] <= y)

    optimize!(data_center_model)

    return value.(p), value.(y)
end

function Google_CAC_No_Peak(emission_rate::Vector{Float64}, emission_cost::Float64, peak_capacity::Float64, flex_load::Vector{Float64}, fixed_load::Vector{Float64}, optimizer)
    data_center_model = Model(optimizer.Optimizer)

    @variable(data_center_model, p[1:length(flex_load)] >= 0)

    @objective(data_center_model, Min, emission_cost * sum(emission_rate[h] * p[h] for h in 1:length(flex_load)))

    @constraint(data_center_model, Eq_total_demand, sum(p[h] for h in 1:length(flex_load)) == sum(flex_load[h] for h in 1:length(flex_load)))
    @constraint(data_center_model, Eq_peak_demand[h in 1:length(flex_load)], p[h] + fixed_load[h] <= peak_capacity)

    optimize!(data_center_model)

    return value.(p)
end