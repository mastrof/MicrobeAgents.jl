# Translation of a microbe constrained by obstacle surfaces.

const MAXITER_TRACE = 16

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
    L = cutoff - R # step length
    convex = isconvex(p)
    τ = convex ? ε : NONCONVEX_TOL * ε # contact tolerance (sdf round-off, finite-difference normals)
    if g ≤ τ
        # touching or overlapping. Blocks if the step goes into the surface by more
        # than τ (round-off after a `Slide` projection is not a contact), if
        # overlapping, or, for `stick`, unless clearly moving away.
        n = surface_normal(p, xe)
        into = dot(ŝ, n) * L < (stick ? τ : -τ)
        overlapping = g < -τ
        (into || overlapping) && return g < acc[4] ? (acc[1], k, i, g, n) : acc
        # moving away from a convex surface never hits it; a non-convex one may be
        # met again further along (concave parts): march tolerating a shallow overlap
        convex && return acc
        δ = _minstep(R, L)
        step = _march(p, xe, ŝ, R, L, -δ, τ, δ)
    elseif convex
        # linear bound (the distance along the ray is convex)
        a = dot(ŝ, surface_normal(p, xe))
        step = a < 0 ? g / -a : Inf
    else
        step = _march(p, xe, ŝ, R, L, 0.0, τ, _minstep(R, L))
    end
    return step < acc[1] ? (step, acc[2], acc[3], acc[4], acc[5]) : acc
end

const NONCONVEX_TOL = 1e3
# minimum ray-marching step: non-convex features thinner than this may be missed
_minstep(R, L) = 1e-2 * max(R, L)

# Distance along unit `ŝ` (up to `L`) from `x` at which the gap `sdf(p) - R` first
# enters the band `[h, h + τ]` (`Inf` if it stays above `h + τ`): sphere tracing with
# minimum step `δ`, then bisection if a step overshoots the band.
function _march(p, x, ŝ, R, L, h, τ, δ)
    lo, glo = 0.0, sdf(p, x) - R - h
    while lo < L
        t = min(lo + max(glo, δ), L)
        gt = sdf(p, x + t * ŝ) - R - h
        if gt ≤ τ
            gt ≥ 0 && return t
            hi = t
            for _ in 1:64
                m = (lo + hi) / 2
                gm = sdf(p, x + m * ŝ) - R - h
                0 ≤ gm ≤ τ && return m
                gm > τ ? (lo = m) : (hi = m)
            end
            return lo
        end
        lo, glo = t, gt
    end
    return Inf
end

# contact rule of set `k` (type-stable tuple lookup)
@inline _rule(sets::Tuple, k) = k == 1 ? first(sets).contact : _rule(Base.tail(sets), k - 1)

# Keep a microbe held by an adhesive rule on the surface it slid along: bring it
# back to contact distance and transport its direction, tangent to the surface at
# the last contact (normal `n0`), to the tangent plane there.
function _readhere(sets, k, microbe, x, n0, R, reach, model)
    for _ in 1:2
        d, n, kn = _nearest(sets, Agents.normalize_position(x, model), R + reach)
        kn == k || return x
        x += (R - d) * n
    end
    _, n, = _nearest(sets, Agents.normalize_position(x, model), R + reach)
    microbe.vel = _transport(direction(microbe), n0, n)
    return x
end

# Direction `v` rotated by the minimal rotation carrying unit normal `n0` onto `n`,
# then made tangent to the surface with normal `n`. The rotation preserves the
# norm: projecting and renormalizing instead inflates the components of `v` off the
# osculating plane, so geodesics drift on surfaces that are not umbilic (cylinders, tori).
# Normals turning by more than 90° in one step (sharp features) only get the projection.
function _transport(v, n0, n)
    a = dot(n0, n)
    if a > 0
        w = n0 + n
        v = v - dot(v, w) / (1 + a) * w + 2dot(v, n0) * n
    end
    t = v - dot(v, n) * n
    l = norm(t)
    return l > 0 ? t / l : v
end

# sliding along concave parts of non-convex solids, or running out of trace
# iterations, may leave a shallow overlap: push the microbe out of it
function _remove_overlap(sets, x, R, model)
    for _ in 1:2
        d, n, = _nearest(sets, Agents.normalize_position(x, model), R)
        d < R || return x
        x += (R - d) * n
    end
    return x
end

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
    nheld = zv # normal at its last contact
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
            rule = _rule(sets, k)
            s = SVector{D,Float64}(contact!(rule, microbe, c, model))
            adheres(rule) && (held = k; nheld = n)
        elseif step ≥ L
            x += s
            break
        else
            x += step * ŝ
            s -= step * ŝ
        end
    end
    held > 0 && (x = _readhere(sets, held, microbe, x, nheld, R, reach, model))
    _has_nonconvex(sets) && (x = _remove_overlap(sets, x, R, model))
    move_agent!(microbe, Agents.normalize_position(x, model), model)
    return microbe
end
