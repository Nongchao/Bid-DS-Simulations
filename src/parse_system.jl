function create_rts_gmlc_system(
    repo;
    time_series_resolution = Dates.Hour(1),
    try_deserialize = true,
    kwargs...,
)

    rts_data = joinpath(repo, "RTS_Data")
    src_data = joinpath(rts_data, "SourceData")
    siip_data = joinpath(rts_data, "FormattedData", "SIIP")
    data = PowerSystemTableData(
        src_data,
        100.0,
        joinpath(siip_data, "user_descriptors.yaml"),
        generator_mapping_file = joinpath(siip_data, "generator_mapping.yaml"),
        timeseries_metadata_file = joinpath(siip_data, "timeseries_pointers.json"),
    )
    sys = System(data; 
        time_series_resolution = time_series_resolution, 
        time_series_directory="/tmp/scratch", 
        kwargs...
    )
    set_units_base_system!(sys, "system_base")

    return sys
end


function add_emission_cost(emission_cost_df, source_gen_data)
    for r in 1:size(source_gen_data)[1]
        for e in emission_type
            if source_gen_data[r, "Emissions " * e * " Lbs/MMBTU"] != "Unit-specific"
                if typeof(source_gen_data[r, "Emissions " * e * " Lbs/MMBTU"]) == String
                    emission_rate = parse(Float64, source_gen_data[r, "Emissions " * e * " Lbs/MMBTU"])
                elseif (typeof(source_gen_data[r, "Emissions " * e * " Lbs/MMBTU"]) == Int64) || (typeof(source_gen_data[r, "Emissions " * e * " Lbs/MMBTU"]) == Float64)
                    emission_rate = source_gen_data[r, "Emissions " * e * " Lbs/MMBTU"]
                else
                    emission_rate = 0.0
                end
            else
                emission_rate = 0.0
            end
            source_gen_data[r, "Fuel Price \$/MMBTU"] = source_gen_data[r, "Fuel Price \$/MMBTU"] + emission_rate * emission_cost_df[1, Symbol(e)]
        end
    end

    CSV.write(joinpath(repo, "RTS_Data", "SourceData", "gen.csv"), source_gen_data)
end
