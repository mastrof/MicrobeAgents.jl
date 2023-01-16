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
function unfold!(unfolded, cnf₁, cnf₀, L)
    dim = length(first(cnf₁))
    nmicrobes, = size(cnf₁)
    for i in 1:nmicrobes
        l = typeof(L) <: AbstractArray ? L[μ] : L
        newx = ntuple(μ -> unfold_coord(cnf₀[i][μ], cnf₁[i][μ], l), dim)
        unfolded[i] = newx
    end # for
end # function

"""
    unfold(trajectory::T, L) where {S<:Tuple, T<:AbstractArray{S,2}}
Unfold `trajectory` in a domain of periodicity `L`.
"""
function unfold(trajectory::AbstractMatrix{<:Tuple}, L)
    nmicrobes, nsteps = size(trajectory)
    unfolded = Matrix{eltype(trajectory)}(undef, size(trajectory)...)
    unfolded[:,1] .= trajectory[:,1]
    for t in 2:nsteps
        oldcnf = unfolded[:,t-1]
        newcnf = trajectory[:,t]
        unfolded_slice = @view unfolded[:,t]
        unfold!(unfolded_slice, newcnf, oldcnf, L)
    end # for
    return unfolded
end # function

"""
    msd(adf; tstep::Int=1, L=Inf)
Evaluate mean-squared displacement from an agent dataframe `adf` containing
the position timeseries of agents (`:pos` field).
Assumes that sampled positions are uniformly spaced in time.
Parameter `L` defines the periodicity of the domain for unfolding;
set `L=Inf` (default) if boundaries are not periodic.
    msd(trajectory::T, tstep::Int=1) where {S,T<:AbstractArray{S,2}}
Evaluate mean-squared displacement of `trajectory`, where different microbes
are collected along first dimension, and times along second dimension.
"""
function msd(adf; tstep::Int=1, L=Inf)
    trajectory = vectorize_adf_measurement(adf, :pos)
    if isinf(L)
        return msd(trajectory, tstep)
    else
        trajectory_unfolded = unfold(trajectory, L)
        return msd(trajectory_unfolded, tstep)
    end # if
end # function

function msd(trajectory::AbstractMatrix{<:Tuple}, tstep::Int=1)
    nmicrobes, nsteps = size(trajectory)
    timelags = range(1, nsteps-1; step=tstep)
    MSD = zeros(nsteps-1)
    for (j,Δt) in enumerate(timelags)
        for t₀ in 1:nsteps-Δt
            for i in 1:nmicrobes
                u = trajectory[i, t₀]
                v = trajectory[i, t₀+Δt]
                MSD[j] += sum(abs2.(u .- v))
            end # for
        end # for
        MSD[j] /= (nmicrobes * (nsteps-j))
    end # for
    return MSD
end # function