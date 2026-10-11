export SurfaceContact, contact!, Halt, Slide, Reflect, Randomize

"""
    SurfaceContact{D}
Contact event passed to `contact!`. Fields:
- `set::Int`: index of the obstacle set in `model.obstacles`
- `index::Int`: index of the primitive within the set
- `point::SVector{D,Float64}`: contact point on the surface
- `normal::SVector{D,Float64}`: outward unit normal of the surface
- `displacement::SVector{D,Float64}`: displacement still to perform in this step
"""
struct SurfaceContact{D}
    set::Int
    index::Int
    point::SVector{D,Float64}
    normal::SVector{D,Float64}
    displacement::SVector{D,Float64}
end

"""
    contact!(rule, microbe, contact::SurfaceContact, model) -> SVector
Apply the contact `rule` of an obstacle set when `microbe` touches its surface.
Returns the displacement still to perform in the current step;
the rule may also modify the microbe (direction, speed, motile state, ...).
A plain function used as rule is called as `rule(microbe, contact, model)`.
"""
contact!(f::Function, microbe, c::SurfaceContact, model) = f(microbe, c, model)

"""
    Halt()
Contact rule: the microbe stops at the surface for the rest of the step,
keeping its direction (it leaves once its direction points away from the surface).
"""
struct Halt end
contact!(::Halt, microbe, c::SurfaceContact, model) = zero(c.displacement)

"""
    Slide(; align = false)
Contact rule: the microbe moves along the surface (normal component of the
displacement removed). With `align = true` the microbe stays on the surface
and keeps its direction tangent to it, following curved surfaces, until
it is reoriented away from the surface.
"""
Base.@kwdef struct Slide
    align::Bool = false
end
# adhesive rules keep the microbe on the surface across steps (see `constrained_move!`)
adheres(rule) = false
adheres(rule::Slide) = rule.align
function contact!(rule::Slide, microbe, c::SurfaceContact, model)
    s, n = c.displacement, c.normal
    a = dot(s, n)
    t = a < 0 ? s - a * n : s
    if rule.align
        l = norm(t)
        l > 0 && (microbe.vel = t / l)
    end
    return t
end

"""
    Reflect()
Contact rule: specular reflection of direction and displacement on the surface.
"""
struct Reflect end
_reflect(v, n) = (a = dot(v, n); a < 0 ? v - 2a * n : v)
function contact!(::Reflect, microbe, c::SurfaceContact, model)
    microbe.vel = _reflect(direction(microbe), c.normal)
    return _reflect(c.displacement, c.normal)
end

"""
    Randomize()
Contact rule: the microbe takes a random direction, uniformly distributed
away from the surface, and continues along it.
"""
struct Randomize end
function contact!(::Randomize, microbe::AbstractMicrobe{D}, c::SurfaceContact, model) where {D}
    u = _reflect(random_velocity(abmrng(model), D), c.normal)
    microbe.vel = u
    return norm(c.displacement) * u
end
