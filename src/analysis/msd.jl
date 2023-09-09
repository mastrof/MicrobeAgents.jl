export msd

"""
    unfold_coord(x₀, x₁, L)
Unfold a translation `x₀` ↦ `x₁` in a domain of periodicity `L`.
"""
function unfold_coord(x₀, x₁, L)
    dx = x₁ - x₀
    sdx = sign(dx)
    a = round(abs(dx/L))
    if abs(dx) > L/2
        return x₁ - a*L*sdx
    else
        return x₁
    end # if
end # function

"""
    unfold!(unfolded, cnf₁, cnf₀, L)
Unfold spatial configuration `cnf₁` with respect to `cnf₀` in a domain of
periodicity `L` and store to `unfolded`.
"""
function unfold!(unfolded, cnf₁, cnf₀, L::Real)
    dim = length(first(cnf₁))
    nmicrobes, = size(cnf₁)
    for i in 1:nmicrobes
        newx = ntuple(μ -> unfold_coord(cnf₀[i][μ], cnf₁[i][μ], L), dim)
        unfolded[i] = newx
    end # for
end # function
function unfold!(unfolded, cnf₁, cnf₀, L::NTuple)
    dim = length(first(cnf₁))
    nmicrobes, = size(cnf₁)
    for i in 1:nmicrobes
        newx = ntuple(μ -> unfold_coord(cnf₀[i][μ], cnf₁[i][μ], L[μ]), dim)
        unfolded[i] = newx
    end 
end

"""
    unfold(trajectory::T, L) where {S<:Tuple, T<:AbstractArray{S,2}}
Unfold `trajectory` in a domain of periodicity `L`.
"""
function unfold(trajectory::AbstractMatrix{<:SVector}, L)
    nsteps, nmicrobes = size(trajectory)
    unfolded = Matrix{eltype(trajectory)}(undef, size(trajectory)...)
    unfolded[1,:] .= trajectory[1,:]
    for t in 2:nsteps
        oldcnf = unfolded[t-1,:]
        newcnf = trajectory[t,:]
        unfolded_slice = @view unfolded[t,:]
        unfold!(unfolded_slice, newcnf, oldcnf, L)
    end # for
    return unfolded
end # function

"""
    msd(adf, lags=1:last(adf.step)-1; L=Inf)
Evaluate mean-squared displacement from an agent dataframe `adf` containing
the position timeseries of agents (`adf.pos`).
Parameter `L` defines the periodicity of the domain for unfolding;
set `L=Inf` (default) if boundaries are not periodic.
"""
function msd(adf, lags=1:last(adf.step)-1; L=Inf)
    trajectory = vectorize_adf_measurement(adf, :pos)
    if L isa Real && isinf(L)
        return msd(trajectory, lags)
    else
        trajectory_unfolded = unfold(trajectory, L)
        return msd(trajectory_unfolded, lags)
    end # if
end # function

"""
    msd(x, lags=1:size(x,1)-1)
Evaluate the mean squared displacement of `x` at timelags `lags`.
If `lags` is unspecified, all available timepoints are used (lag-0 excluded).

`x` can be an `AbstractVector` or an `AbstractMatrix`. In the latter case,
it is assumed that each column represents a distinct signal and each row
a timepoint. The msd is then evaluated for each individual timeseries and
the results are averaged together.
"""
function msd(x::AbstractMatrix, lags=1:size(x,1)-1)
    mean(mapslices(y -> msd(y, lags), x, dims=1), dims=2)
end
function msd(x::AbstractVector, lags=1:size(x,1)-1)
    l = size(x,1)
    S₂ = acf(x,lags)
    D = OffsetArray([0; dot.(x,x); 0], -1:l)
    Q = 2*sum(D)
    S₁ = AxisArray(similar(S₂), lags)
    for k in 1:lags[end]
        Q -= (D[k-1] + D[l-k])
        if k ∈ lags
            S₁[atvalue(k)] = Q / (l-k)
        end
    end
    @. S₁-2S₂
end
