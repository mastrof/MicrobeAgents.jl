function turn!(microbe::AbstractMicrobe{1}, model)
    microbe.vel = rotate(direction(microbe))
end
function turn!(microbe::AbstractMicrobe, model)
    e = direction(microbe)
    θ = rand(abmrng(model), angle(motilepattern(microbe)))
    φ = rand(abmrng(model), azimuthal(motilepattern(microbe)))
    microbe.vel = rotate(e, θ, φ)
end

"""
    rotational_diffusion!(microbe, model)
Reorient `microbe` due to brownian rotational diffusion.
In 1-dimensional models, this functions does nothing.
"""
rotational_diffusion!(microbe::AbstractMicrobe{1}, model) = nothing
function rotational_diffusion!(microbe::AbstractMicrobe{2}, model)
    dt = model.timestep
    D_rot = rotational_diffusivity(microbe)
    σ = sqrt(2*D_rot*dt)
    θ = rand(abmrng(model), Normal(0, σ))
    microbe.vel = rotate(direction(microbe), θ)
end
function rotational_diffusion!(microbe::AbstractMicrobe{3}, model)
    dt = model.timestep
    D_rot = rotational_diffusivity(microbe)
    σ = sqrt(2*D_rot*dt)
    θ = rand(abmrng(model), Normal(0,σ))
    φ = rand(abmrng(model), Uniform(-π, +π))
    microbe.vel = rotate(direction(microbe), θ, ϕ)
    nothing
end

"""
    rotate(w::SVector{D}, θ, ϕ) where D
Rotate `D`-dimensional vector `w` by angles `θ` and `ϕ`.
`θ` is the polar (or inclination) angle, `ϕ` is the azimuthal angle
around the original direction of `w`.

Supports 1, 2, and 3 dimensions.
In `D=1` the rotation is just an inversion, in `D=2` the rotation is only
described by `θ`.
"""
rotate(w::SVector{1}, θ, ϕ) = rotate(w)
rotate(w::SVector{1}, θ) = rotate(w)
rotate(w::SVector{1}) = -w

rotate(w::SVector{2}, θ, ϕ) = rotate(w, θ)
function rotate(w::SVector{2}, θ)
    s, c = sincos(θ)
    SVector{2}(c*w[1]-s*w[2], s*w[1]+c*w[2])
end

function rotate(w::SVector{3}, θ, ϕ)
    # find one vector u which is normal to w
    # rotate w around its normal u by θ
    # then rotate the resulting vector around w by ϕ
    u = normalvector(w)
    q0 = Quaternion(0, w[1], w[2], w[3])
    q1 = quaternion_from_axisangle(u, θ)
    q2 = quaternion_from_axisangle(w, ϕ)
    q = (q2*q1) * q0 * conj(q2*q1)
    SVector{3}(imag_part(q))
end
"""
    normalvector(w::SVector{3})
Find a vector normal to `w`.
"""
function normalvector(w::SVector{3})
    m = findfirst(w .≠ 0)
    n = m%3 + 1
    u = SVector(0., 0., 0.)
    u = setindex(u, w[m], n)
    u = setindex(u, -w[n], m)
end

"""
    quaternion_from_axisangle(axis, angle)
Obtain a quaternion representation of the rotation around `axis` by `angle`.
"""
function quaternion_from_axisangle(axis::SVector{3,<:Real}, angle::Real)
    s, c = sincos(angle / 2)
    axis = normalize(axis)
    return Quaternion(c, s*axis[1], s*axis[2], s*axis[3])
end
