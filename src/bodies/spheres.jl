export HyperSphere, contact, is_encounter

# dispatch hides call to Point
HyperSphere(center::SVector{D}, radius::Real) where D = HyperSphere(Point(Float64.(center)), Float64(radius))

@inline Base.position(a::HyperSphere{D}) where D = SVector{D}(a.center)
#@inline radius(a::AbstractMicrobe) = a.radius
@inline radius(a::HyperSphere) = a.r
@inline contact(a,b,model) = distance(a,b,model) ≤ radius(a) + radius(b)

"""
    is_encounter(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
Check if an encounter occurred between `microbe` and `sphere`.
"""
function is_encounter(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
    d = distance(microbe, sphere, model)
    R = radius(microbe) + radius(sphere)
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
    R = radius(microbe) + radius(sphere)
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
