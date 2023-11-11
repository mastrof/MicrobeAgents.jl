export imsd!, imsd, emsd

function imsd!(df::GroupedDataFrame, sym)
    transform!(df, sym => imsd)
end
function MeanSquaredDisplacement.imsd(df::AbstractDataFrame, sym)
    imsd(groupby(df, :id), sym)
end
function MeanSquaredDisplacement.imsd(df::GroupedDataFrame, sym)
    imsd.([g[!,sym] for g in df])
end

function MeanSquaredDisplacement.emsd(df::AbstractDataFrame, sym)
    emsd(groupby(df, :id), sym)
end
function MeanSquaredDisplacement.emsd(df::GroupedDataFrame, sym)
    emsd([g[!,sym] for g in df])
end
