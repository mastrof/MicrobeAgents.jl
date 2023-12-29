function turn!(microbe::AbstractMicrobe, model)
    e = direction(microbe)
    # extract new speed and rotation angles
    U₁, θ, ϕ = rand(abmrng(model), motilepattern(microbe))
    # reorient and update direction
    microbe.vel = rotate(e, θ, ϕ)
    microbe.speed = U₁
    # switch motile state (does nothing if motility is one-step)
    switch!(motilepattern(microbe))
end

"""
    rotational_diffusion!(microbe, model)
Reorient `microbe` due to brownian rotational diffusion.
In 1-dimensional models, this functions does nothing.
"""
rotational_diffusion!(microbe::AbstractMicrobe{1}, model) = nothing
function rotational_diffusion!(microbe::AbstractMicrobe, model)
    dt = model.timestep
    D_rot = rotational_diffusivity(microbe)
    σ = sqrt(2*D_rot*dt)
    θ = rand(abmrng(model), Normal(0,σ))
    ϕ = rand(abmrng(model), Arccos())
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
rotate(w::SVector{2}, θ) = Angle2d(θ)*w

function rotate(w::SVector{3}, θ, ϕ)
    # find one vector u which is normal to w
    # rotate w around its normal u by θ
    # then rotate the resulting vector around w by ϕ
    u = normalvector(w)
    AngleAxis(ϕ, w...) * AngleAxis(θ, u...) * w
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
