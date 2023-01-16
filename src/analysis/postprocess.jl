export vectorize_adf_measurement

"""
    vectorize_adf_measurement(adf, sym)
Collect quantity `sym` from the agent dataframe `adf` and return it in matrix 
form with dimensions (microbes, times).
"""
function vectorize_adf_measurement(adf, sym)
    nmicrobes = unique(adf[!,:id]) |> length
    nsteps = unique(adf[!,:step]) |> length
    datatype = typeof(adf[1,sym])
    s = Matrix{datatype}(undef, nmicrobes, nsteps)
    for t in 1:nsteps
        for i in 1:nmicrobes
            s[i,t] = adf[i + (t-1)*nmicrobes, sym]
        end # for
    end # for
    return s
end # function