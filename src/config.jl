#using ConfParser
#using Distributions
#include("spot.jl")
#include("evolution.jl")

DEG2RAD = pi / 180.0

function readin(filename::String)
    conf = ConfParse(filename)
    parse_conf!(conf)
    return conf
end

function conf_spotmodel(filename::String)
    conf = readin(filename)
    nspot = retrieve(conf, "spots", "nspot", Int)
    rstar = retrieve(conf, "star", "radius", Float64)
    period = ones(Number, nspot) .* retrieve(conf, "star", "period", Float64)
    u = ones(Number, nspot) .* retrieve(conf, "star", "u", Float64)
    cspot = ones(Number, nspot) .* retrieve(conf, "spots", "cspot", Float64)
    vconv = retrieve(conf, "star", "vconv", Float64)
    cfac = ones(Number, nspot) .* retrieve(conf, "spots", "cfac", Float64)
    amax = ones(Number, nspot) .* retrieve(conf, "star", "amax", Float64)
    pk = ones(Number, nspot)
    phase = rand(Uniform(0.0, 2 * pi), nspot) # initial random
    incl = ones(Number, nspot) .* (retrieve(conf, "star", "incl", Float64) * DEG2RAD)
    if haskey(conf._data["spots"], "lat")#set lat according to pk
        lat = retrieve(conf, "spots", "lat", Float64)
        lat = rand(Uniform(-lat, lat), nspot) .* DEG2RAD
    else
        lat = zeros(Number, nspot) .* DEG2RAD
    end
    Q = ones(Number, nspot) .* retrieve(conf, "spots", "Q", Float64)

    if retrieve(conf, "spots", "decay", Bool)
        decay = period .* 1.0 .^ rand(Normal(0, 0.5), nspot)
        #need complete the snippet to read in decay block, here we can consider multiple dispatch
    else
        decay = zeros(Number, nspot)
    end

    if haskey(conf, "diagram")
        t0 = 0
        cycle = retrieve(conf, "diagram", "cycle", Float64)
        init_lat = retrieve(conf, "diagram", "init_lat", Float64)
        diagram = Diagram(t0, cycle, init_lat)
        SpotModel = RotationSpot(nspot, rstar, period, u, cspot, vconv, cfac,
        amax, pk, decay, phase, incl, lat, Q)
        if retrieve(conf, "spots", "differential", Float64) != 0
            diffrot = retrieve(conf, "spots", "differential", Float64)
            omega_eq = 2 * pi / retrieve(conf, "star", "period", Float64)
            differential = Differential(omega_eq, diffrot)

            return [SpotModel, diagram, differential]
        else
            return [SpotModel, diagram]
        end
    else
        SpotModel = RotationSpot(nspot, rstar, period, u, cspot, vconv, cfac,
        amax, pk, decay, phase, incl, lat, Q)
        if retrieve(conf, "spots", "differential", Float64) != 0
            diffrot = retrieve(conf, "spots", "differential", Float64)
            omega_eq = 2 * pi / retrieve(conf, "star", "period", Float64)
            differential = Differential(omega_eq, diffrot)

            return [SpotModel, differential]
        else
            return SpotModel
        end
    end
end
