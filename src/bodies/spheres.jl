export HyperSphere, contact, is_encounter

# dispatch hides call to Point
HyperSphere(center::NTuple{D}, radius::Real) where D = HyperSphere(Point(Float64.(center)), Float64(radius))

distance(a::HyperSphere{D}, b::HyperSphere{D}, model) where D = distance(Tuple(a.center), Tuple(b.center), model)
distance(a::AbstractMicrobe{D}, b::HyperSphere{D}, model) where D = distance(a.pos, Tuple(b.center), model)
distance(a::HyperSphere{D}, b::AbstractMicrobe{D}, model) where D = distance(Tuple(a.center), b.pos, model)
distance(a::NTuple{D}, b::HyperSphere{D}, model) where D = distance(a, Tuple(b.center), model)
distance(a::HyperSphere{D}, b::NTuple{D}, model) where D = distance(Tuple(a.center), b, model)

contact(a::HyperSphere{D}, b::AbstractMicrobe{D}, model) where D = distance(a,b,model) ≤ a.r+b.radius
contact(a::AbstractMicrobe{D}, b::HyperSphere{D}, model) where D = distance(a,b,model) ≤ a.radius+b.r
contact(a::HyperSphere{D}, b::HyperSphere{D}, model) where D = distance(a,b,model) ≤ a.r+b.r

"""
    is_encounter(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
Check if an encounter is occurring between `microbe` and `sphere`.
"""
function is_encounter(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
    # if microbe and sphere are in contact, return true immediately
    if contact(microbe, sphere, model)
        return true
    else # check if an intersection will occur during the next step
        return line_sphere_intersection(microbe, sphere, model)
    end
end

"""
    line_sphere_intersection(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
Check whether `microbe` will intersect `sphere` during its next step.
"""
function line_sphere_intersection(microbe::AbstractMicrobe{D}, sphere::HyperSphere{D}, model)::Bool where D
    c = sphere.center
    x1 = @. microbe.pos - c
    x2 = @. x1 + microbe.vel*microbe.speed*model.timestep
    R = microbe.radius + sphere.r
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
