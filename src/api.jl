export position, direction, speed, velocity, motilepattern,
    turnrate, rotational_diffusivity, radius, state,
    distance, distancevector

"""
    position(m::AbstractMicrobe)
Return the position of the microbe.
"""
Base.position(m::AbstractMicrobe) = m.pos
"""
    direction(m::AbstractMicrobe)
Return the direction versor of the microbe.
"""
direction(m::AbstractMicrobe) = m.vel
"""
    speed(m::AbstractMicrobe)
Return the speed of the microbe.
"""
speed(m::AbstractMicrobe) = m.speed
"""
    velocity(m::AbstractMicrobe)
Return the velocity vector of the microbe (direction times speed).
"""
velocity(m::AbstractMicrobe) = direction(m) .* speed(m)
"""
    motilepattern(m::AbstractMicrobe)
Return the motile pattern of the microbe.
"""
motilepattern(m::AbstractMicrobe) = m.motility
"""
    turnrate(m::AbstractMicrobe)
Return the unbiased turn rate of the microbe.
"""
turnrate(m::AbstractMicrobe) = m.turn_rate
"""
    rotational_diffusivity(m::AbstractMicrobe)
Return the rotational diffusivity of the microbe.
"""
rotational_diffusivity(m::AbstractMicrobe) = m.rotational_diffusivity
"""
    radius(m::AbstractMicrobe)
Return the radius of the microbe.
"""
radius(m::AbstractMicrobe) = m.radius
"""
    state(m::AbstractMicrobe)
Return the internal state of the microbe.
"""
state(m::AbstractMicrobe) = m.state

"""
    distance(a, b, model)
Evaluate the euclidean distance between `a` and `b` respecting
the boundary conditions of the `model`
"""
distance(a, b, model) = euclidean_distance(position(a), position(b), model)
"""
    distancevector(a, b, model)
Evaluate the distance vector from `a` to `b` respecting
the boundary conditions of the `model`.
"""
distancevector(a, b, model) = distancevector(position(a), position(b), model)
function distancevector(a::SVector{D}, b::SVector{D}, model) where D
    extent = spacesize(model)
    SVector{D}(wrapcoord(a[i], b[i], extent[i]) for i in 1:D)
end

Base.position(a::SVector{D}) where D = a
Base.position(a::NTuple{D}) where D = SVector{D}(a)

## utils
function wrapcoord(x1, x2, d)
    a = (x2 - x1) / d
    (a - round(a)) * d
end
