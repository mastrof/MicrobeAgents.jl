export vectorize_adf_measurement

"""
    vectorize_adf_measurement(adf, sym)
Collect quantity `sym` from the agent dataframe `adf` and return it in matrix 
form with dimensions (times, microbes).
"""
function vectorize_adf_measurement(adf, sym)
    nmicrobes = unique(adf[!,:id]) |> length
    nsteps = unique(adf[!,:step]) |> length
    datatype = typeof(adf[1,sym])
    s = Matrix{datatype}(undef, nsteps, nmicrobes)
    for i in 1:nmicrobes
        for t in 1:nsteps
            s[t,i] = adf[i + (t-1)*nmicrobes, sym]
        end # for
    end # for
    return s
end # function