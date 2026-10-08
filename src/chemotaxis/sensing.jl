"""
Conversion factor (1/√(number of molecules) --> 1/√(moles)) used
in the evaluation of chemotactic sensing noise.
"""
global const CONV_NOISE::Float64 = 0.04075

# Berg-Purcell noise scales as 1/√radius, so noisy sensing needs a finite radius
function check_sensing_radius(behavior, Π, microbe)
    if Π > 0 && iszero(radius(microbe))
        throw(ArgumentError(
            "$(nameof(typeof(behavior))) with chemotactic_precision > 0 requires " *
            "a microbe radius > 0 (e.g. `add_agent!(model; radius = 0.5, ...)`)"
        ))
    end
    return nothing
end
