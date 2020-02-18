#using Distributions

mutable struct Diagram
    t0::Float64
    cycle::Float64
    init_lat::Float64
end

struct Differential
    omega_eq::Float64
    diffrot::Float64
end

function sample_lat(diagram::Diagram, t::Float64)
    if (t - diagram.t0) <= 0
        phase = (t - diagram.t0) / (diagram.cycle) % 1
        ave_lat = -1.0 * diagram.init_lat * phase
    else
        phase = (t - diagram.t0) / (diagram.cycle) % 1
        ave_lat = -1.0 * diagram.init_lat * phase + diagram.init_lat
    end
    BiGaussian = MixtureModel(map(u -> Normal(u, 2.5), [ave_lat, -1.0*ave_lat]))

    return rand(BiGaussian)
end

function sample_lifetime(decay_scale::Float64)
    u = MixtureModel(Uniform.([2/24*decay_scale, 2.0*decay_scale, 30.0*decay_scale], 
                            [2.0*decay_scale, 30.0*decay_scale, 60.0*decay_scale]),
                            [0.0, 0.2, 0.8])
    return rand(u)
end

function sample_phase(lat::Float64) # TODO: consider the effect of shear force?
    if lat >= 0
        binor = MixtureModel(map(u -> Normal(u, 15 * pi / 180), [0, pi]))
        return rand(binor)
    else
        binor = MixtureModel(map(u -> Normal(u, 15 * pi / 180), [pi / 2, 3 * pi / 2]))
        return rand(binor)
    end
end
