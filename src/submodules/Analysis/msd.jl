export imsd!, emsd

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
        newx = SVector{dim}(unfold_coord(cnf₀[i][μ], cnf₁[i][μ], L) for μ in 1:dim)
        unfolded[i] = newx
    end # for
end # function
function unfold!(unfolded, cnf₁, cnf₀, L::SVector)
    dim = length(first(cnf₁))
    nmicrobes, = size(cnf₁)
    for i in 1:nmicrobes
        newx = SVector{dim}(unfold_coord(cnf₀[i][μ], cnf₁[i][μ], L[μ]) for μ in 1:dim)
        unfolded[i] = newx
    end 
end

"""
    unfold(trajectory::T, L) where {S<:SVector, T<:AbstractArray{S,2}}
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

function imsd!(df::GroupedDataFrame, sym)
    transform!(df, sym => imsd)
end

function MeanSquaredDisplacement.emsd(df::GroupedDataFrame, sym)
    emsd([g.vel for g in df])
end
