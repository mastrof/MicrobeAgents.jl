export SphericalSurfaces

module SphericalSurfaces

using MicrobeAgents
using Agents
using CellListMap.PeriodicSystems
using LinearAlgebra: dot, norm
using Rotations
using GeometryBasics: HyperSphere, Point

export HyperSphere, contact
export surface_interactions!, microbe_surface_step!
export stick!, slide!, stickyslide!
export is_encounter


#== Utility functions ==#
# dispatch hides call to Point
"""
    HyperSphere(center::NTuple{D,<:Real}, radius::Real) where D
Create a `D`-dimensional sphere with radius `radius` and origin at `center`.
"""
HyperSphere(center::NTuple{D}, radius::Real) where D = HyperSphere(Point(Float64.(center)), Float64(radius))

# define _pos to interface with distance functions
using MicrobeAgents: _pos
@inline MicrobeAgents._pos(a::HyperSphere{D}) where D = Tuple(a.center)

@inline _radius(a::NTuple{D,<:Real}) where D = 0.0
@inline _radius(a::AbstractMicrobe) = a.radius
@inline _radius(a::HyperSphere{D}) where D = a.r
"""
    contact(a, b, model)
Check whether `a` and `b` are at contact, given the spatial properties of model.
`a` and `b` can be any mix of `AbstractMicrobe`s, `HyperSphere`s and `NTuple`s.
For `NTuple`s, it is assumed that the points have radius 0.
"""
@inline contact(a,b,model) = distance(a,b,model) ≲ _radius(a) + _radius(b)

@inline safe_acos(x::T) where T = x ≈ 1 ? zero(x) : x ≈ -1 ? T(π) : acos(x)
@inline ≲(x, y) = x ≤ y || isapprox(x, y; atol=1e-10) # \lesssim


#== Surface Interactions ==#
"""
    surface_interactions!(model, spheres, interaction=stickyslide!)
Setup `model` to evaluate surface interactions between its microbes and `spheres`.

Three types of surface interactions are available out-of-the-box:
`stick!`, `slide!` and `stickyslide!`.
The `interaction` defaults to `stickyslide!`
"""
function surface_interactions!(
    model::ABM, spheres::AbstractVector, interaction=stickyslide!
)
    abmproperties(model)[:spheres] = spheres
    abmproperties(model)[:surface!] = interaction
    #== used only for stickyslide! ==#
    abmproperties(model)[:is_stuck] = fill(false, nagents(model))
    abmproperties(model)[:slidingdirection] = fill(+1, nagents(model))
end


#== these are not used at the moment
function surface_interaction!(model::ABM, listkey::Symbol)
    map_pairwise!(
        (x, y, i, j, d², f) -> surface_interaction!(x, y, i, j, d², f, model),
        abmproperties(model)[listkey]
    )
end
function surface_interaction!(x, y, i, j, d², f, model)
    microbe = model[i]
    sphere = abmproperties(model)[:spheres][j]
    model.surface_affect!(microbe, sphere, model)
end
==#


"""
    microbe_surface_step!(microbe, model)
Perform an integration step for `microbe`, executing the following operations:
1. Update microbe position according to its current velocity
2. Reorient microbe through rotational diffusion
3. Apply surface interactions (as specified by `model.surface!`)
4. Update the state of the microbe through the dedicated `affect!` function
5. Evaluate instantaneous turn rate through `turnrate` function
6. Reorient microbe if necessary

Three types of surface interactions are available out-of-the-box:
`stick!`, `slide!` and `stickyslide!`.
If `stickyslide!` is used, step 1 and 2 are performed only if the microbe is
not stuck to a surface; instead, when the microbe is stuck, it will just
slide smoothly along the sphere surface with constant velocity until
it detaches by performing a tumble.
"""
function microbe_surface_step!(microbe, model)
    dt = model.timestep
    if ~model.is_stuck[microbe.id]
        move_agent!(microbe, model, microbe.speed*dt)
        MicrobeAgents.rotational_diffusion!(microbe, model)
    end
    for sphere in model.spheres
        model.surface!(microbe, sphere, model)
    end
    affect!(microbe, model)
    ω = turnrate(microbe, model)
    if rand(abmrng(model)) < ω*dt
        MicrobeAgents.turn!(microbe, model)
        model.is_stuck[microbe.id] = false
    end
end


"""
    stick!(microbe, sphere, model)
`microbe` sticks to the surface of `sphere` at the point of contact.
"""
function stick!(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model) where D
    R = _radius(microbe) + _radius(sphere)
    R² = R*R
    d = distancevector(microbe, sphere, model)
    d² = sum(abs2.(d))
    c = d² - R²
    if c < 0
        # step the microbe backward to the surface of the sphere
        Δt = model.timestep
        s = @. -microbe.vel*microbe.speed*Δt
        a = sum(abs2.(s))
        b = -2*dot(d, s)
        ε = - b/2a * (1 - sqrt(1 - 4a*c/b^2))
        z = ε .* s
        walk!(microbe, z, model)
    end
end

"""
    slide!(microbe, sphere, model)
`microbe` slides smoothly along the surface of `sphere`.
"""
function slide!(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model) where D
    R = _radius(microbe) + _radius(sphere)
    d = distancevector(microbe, sphere, model)
    Δ = sqrt(sum(abs2.(d)))
    if Δ < R
        ε = 1 - R/Δ
        z = ε .* d
        walk!(microbe, z, model)
    end
end

"""
    stickyslide!(microbe, sphere, model)
`microbe` keeps sliding along the surface of `sphere` until a tumble occurs.
"""
function stickyslide!(microbe::AbstractMicrobe{2}, sphere::HyperSphere{2}, model)
    if ~model.is_stuck[microbe.id]
        stick!(microbe, sphere, model)
        model.is_stuck[microbe.id] = contact(microbe, sphere, model)
        model.slidingdirection[microbe.id] = slidedirection(microbe, sphere, model)
        return microbe # type stability
    end
    s = if ~haskey(abmproperties(model), :slidemode)
        model.slidingdirection[microbe.id]
    elseif uppercase(model.slidemode) == "CCW"
        +1
    elseif uppercase(model.slidemode) == "CW"
        -1
    end
    Δt = model.timestep
    R = _radius(microbe) + _radius(sphere)
    ω = s * microbe.speed / R # angular velocity
    d = distancevector(sphere, microbe, model)
    θ = atan(d[2], d[1]) + ω*Δt
    δ = (cos(θ), sin(θ))
    new_pos = _pos(sphere) .+ R.*δ
    move_agent!(microbe, new_pos, model)
end
function stickyslide!(microbe::AbstractMicrobe{3}, sphere::HyperSphere{3}, model)
    if ~model.is_stuck[microbe.id]
        stick!(microbe, sphere, model)
        model.is_stuck[microbe.id] = true
        model.slidingdirection[microbe.id] = slidedirection(microbe, sphere, model)
        return microbe # type stability
    else
        s = if ~haskey(abmproperties(model), :slidemode)
            model.slidingdirection[microbe.id]
        elseif uppercase(model.slidemode) == "CCW"
            +1
        elseif uppercase(model.slidemode) == "CW"
            -1
        end
        Δt = model.timestep
        R = _radius(microbe) + _radius(sphere)
        ω = s * microbe.speed / sphere.r # angular velocity
        d = distancevector(sphere, microbe, model)
        x, y, z = d
        θ = safe_acos(z/R) + ω*Δt
        ϕ = sign(y)*safe_acos(x / (x^2+y^2)) + ω*Δt
        δ = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
        new_pos = _pos(sphere) .+ R.*δ
        move_agent!(microbe, new_pos, model)
    end
    nothing
end

function slidedirection(microbe::AbstractMicrobe{2}, sphere::HyperSphere{2}, model)
    # rotation_between only works with 3d vectors
    d = [distancevector(sphere, microbe, model)..., 0.0]
    M = rotation_between(d, [1, 0, 0])
    Vy = M[2,1]*microbe.vel[1] + M[2,2]*microbe.vel[2] # y component of rotated vel
    Vy ≥ 0 ? +1 : -1
end


#== Encounters ==#
"""
    is_encounter(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D                 end # module
Check if an encounter occurred between `microbe` and `sphere`.
"""
function is_encounter(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
    d = distance(microbe, sphere, model)
    R = _radius(microbe) + _radius(sphere)
    # if microbe and sphere are in contact, return true immediately
    if d < R
        return true
    # if they are not in contact but close, check if an intersection occurred
    # during the last integration step
    elseif d < R + 2*microbe.speed*model.timestep
        return line_sphere_intersection(microbe, sphere, model)
    else
        return false
    end
end

"""
    line_sphere_intersection(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
Check if `microbe` has intersected `sphere` during the last integration step.
"""
function line_sphere_intersection(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
    x1 = distancevector(model.old_positions[microbe.id], sphere, model)
    x2 = distancevector(microbe, sphere, model)
    R = _radius(microbe) + _radius(sphere)
    line_sphere_intersection(x1,x2,R)
end
"""
    line_sphere_intersection(x1, x2, R)
Check whether a segment between points `x1` and `x2` intersects a sphere
of radius `R` centered at the origin.
"""
function line_sphere_intersection(x1, x2, R)
    dx = x2.-x1
    a = dot(dx,dx)
    b = 2 * dot(dx,x1)
    c = dot(x1,x1) - R*R
    S = b*b - 4*a*c
    if S < 0
        return false
    else
        u1 = (-b + √S) / (2a)
        u2 = (-b - √S) / (2a)
        return (0 ≤ u1 ≤ 1) || (0 ≤ u2 ≤ 1)
    end
end

end # module
