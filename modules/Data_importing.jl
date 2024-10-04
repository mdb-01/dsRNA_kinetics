
"""
    importlabdata()

Import excel file with data from standard at different concentrations

"""
function importlabdata(filename::String)
    # load excel file into DataFrame
    raw_data = CSV.File(filename) |> DataFrame
    raw_data_matrix = Matrix(raw_data)
    
    raw_data_matrix = map(x -> tryparse(Float64, string(x)) !== nothing ? parse(Float64, string(x)) : NaN, raw_data_matrix)

    # create a table with condition number, replicate number, time (hr), mass (kg), NTP conc in Mx4
    processed_data = zeros(size(raw_data_matrix, 1), 2) 

    # fill processed data matrix
    processed_data[:,1] = raw_data_matrix[:,1] # concentration
    processed_data[:,2] = raw_data_matrix[:,2] # intensity (local bc corr volume)

    standard_data_df = DataFrame(processed_data, [:concentration, :volume])

    standard_data_df = standard_data_df[.!isnan.(standard_data_df.concentration) .& .!isnan.(standard_data_df.volume), :]
    
    return standard_data_df
end