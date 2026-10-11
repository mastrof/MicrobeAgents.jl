# Translation of a microbe constrained by obstacle surfaces.

const MAXITER_TRACE = 8

# Scan the candidates of all sets around `x` for a step along unit `ŝ`.
# Returns `(step, k, i, g, n)`: `step` is a safe advance along `ŝ`;
# if `k > 0`, primitive `i` of set `k` blocks the motion (touching with gap `g`
# and approaching, or overlapping), with outward normal `n`.
# `fresh` is true on the first iteration of a step: a microbe held on the surface
# by an adhesive rule (see `adheres`) is then in contact even when moving tangentially.
@inline _scan(::Tuple{}, k, x, ŝ, R, cutoff, ε, fresh, acc) = acc
@inline function _scan(sets::Tuple, k, x, ŝ, R, cutoff, ε, fresh, acc)
    set = first(sets)
    prims = set.primitives
    stick = fresh && adheres(set.contact)
    # closure locals must not shadow outer locals (would box them and allocate)
    op = (a, j, xe) -> _visit(a, prims[j], k, j, xe, ŝ, R, cutoff, ε, stick)
    _scan(Base.tail(sets), k + 1, x, ŝ, R, cutoff, ε, fresh, _reduce_candidates(op, acc, set, x, cutoff))
end

@inline function _visit(acc, p, k, i, xe, ŝ, R, cutoff, ε, stick)
    d = sdf(p, xe)
    d > cutoff && return acc
    g = d - R
    n = surface_normal(p, xe)
    a = dot(ŝ, n)
    if g ≤ ε
        # touching: blocks if the step would move into the surface by more than ε
        # (round-off after a `Slide` projection is not a contact) or if
        # overlapping (or, if `stick`, also when not clearly moving away);
        # otherwise ignored, since moving away from a convex surface never
        # hits it. `cutoff - R` is the step length.
        if (a * (cutoff - R) < (stick ? ε : -ε) || g < -ε) && g < acc[4]
            return (acc[1], k, i, g, n)
        end
        return acc
    end
    # safe step: linear bound for convex primitives (the distance along the ray
    # is convex), sphere tracing otherwise
    step = isconvex(p) ? (a < 0 ? g / -a : Inf) : g
    return step < acc[1] ? (step, acc[2], acc[3], acc[4], acc[5]) : acc
end

# apply the contact rule of set `k`
@inline _contact!(::Tuple{}, k, microbe, c, model) = error("invalid obstacle set index")
@inline _contact!(sets::Tuple, k, microbe, c, model) = k == 1 ?
    contact!(first(sets).contact, microbe, c, model) :
    _contact!(Base.tail(sets), k - 1, microbe, c, model)

# Keep a microbe held by an adhesive rule on the surface it slid along: bring it
# back to contact distance and align its direction with the surface there,
# keeping the sense of travel.
function _readhere(sets, k, microbe, x, R, reach, model)
    for _ in 1:2
        d, n, kn = _nearest(sets, Agents.normalize_position(x, model), R + reach)
        kn == k || return x
        x += (R - d) * n
    end
    _, n, = _nearest(sets, Agents.normalize_position(x, model), R + reach)
    v = direction(microbe)
    t = v - dot(v, n) * n
    l = norm(t)
    l > 0 && (microbe.vel = t / l)
    return x
end

@inline _adheres(::Tuple{}, k) = false
@inline _adheres(sets::Tuple, k) = k == 1 ? adheres(first(sets).contact) : _adheres(Base.tail(sets), k - 1)

_has_nonconvex(sets::Tuple) = any(set -> !isconvex(eltype(set.primitives)), sets)

"""
    constrained_move!(microbe, obstacles::Tuple, model, dt)
Translate `microbe` by `velocity(microbe) * dt` without crossing obstacle surfaces;
on contact, the rule of the touched obstacle set decides how motion continues.
"""
function constrained_move!(microbe::AbstractMicrobe{D}, sets::Tuple, model, dt) where {D}
    abmspace(model).update_vel!(microbe, model)
    R = Float64(radius(microbe))
    ε = 1e-9 * max(R, 1.0)
    x = SVector{D,Float64}(position(microbe))
    s = SVector{D,Float64}(velocity(microbe) .* dt)
    zv = zero(SVector{D,Float64})
    held = 0 # set whose adhesive rule applied during this step
    reach = norm(s)
    for iter in 1:MAXITER_TRACE
        L = norm(s)
        L > 0 || break
        ŝ = s / L
        x = Agents.normalize_position(x, model)
        step, k, i, g, n = _scan(sets, 1, x, ŝ, R, R + L, ε, iter == 1, (Inf, 0, 0, Inf, zv))
        if k > 0
            g < 0 && (x -= g * n) # remove overlap
            c = SurfaceContact{D}(k, i, x - R * n, n, s)
            s = SVector{D,Float64}(_contact!(sets, k, microbe, c, model))
            _adheres(sets, k) && (held = k)
        elseif step ≥ L
            x += s
            break
        else
            x += step * ŝ
            s -= step * ŝ
        end
    end
    held > 0 && (x = _readhere(sets, held, microbe, x, R, reach, model))
    if _has_nonconvex(sets) # motion near non-convex solids is not exact: remove overlaps
        d, n, = _nearest(sets, Agents.normalize_position(x, model), R)
        d < R && (x += (R - d) * n)
    end
    move_agent!(microbe, Agents.normalize_position(x, model), model)
    return microbe
end
