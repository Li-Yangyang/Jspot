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

function corrected_lat(diagram::Diagram, t::Float64)
    if (t - diagram.t0) <= 0
        phase = (t - diagram.t0) / (diagram.cycle * 365.0) % 1
        ave_lat = -1.0 * diagram.init_lat * phase
    else
        phase = (t - diagram.t0) / (diagram.cycle * 365.0) % 1
        ave_lat = -1.0 * diagram.init_lat * phase + diagram.init_lat
    end
    BiGaussian = MixtureModel(map(u -> Normal(u, 2.5), [ave_lat, -1.0*ave_lat]))

    return rand(BiGaussian)
end
