export position, direction, speed, velocity, motilepattern,
    turnrate, rotational_diffusivity, radius, state,
    cwbias,
    distance, distancevector

Base.position(m::AbstractMicrobe) = m.pos
direction(m::AbstractMicrobe) = m.vel
speed(m::AbstractMicrobe) = m.speed
velocity(m::AbstractMicrobe) = direction(m) .* speed(m)
motilepattern(m::AbstractMicrobe) = m.motility
turnrate(m::AbstractMicrobe) = m.turn_rate
rotational_diffusivity(m::AbstractMicrobe) = m.rotational_diffusivity
radius(m::AbstractMicrobe) = m.radius
state(m::AbstractMicrobe) = m.state

distance(a, b, model) = euclidean_distance(position(a), position(b), model)
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
