#using ConfParser
#using Distributions
#include("spot.jl")
include("evolution.jl")

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
        #multinormal distribution according to Dumusque 2016
    else
        decay = zeros(Number, nspot)
    end

    amax = ones(Number, nspot) .* retrieve(conf, "star", "amax", Float64)

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

function initial_condition(Nspot::Int, step::Float64, decay_scale::Float64, t::Array{Float64,1}, diagram::Diagram)
    N = zeros(0)
    decay = zeros(0)
    lat = zeros(0)
    phase = zeros(0)
    pk = zeros(0)
    Ngroup = zeros(0)
    while sum(Ngroup) <= Nspot
        lifetime = sample_lifetime(decay_scale)
        amax = 1e-4 * lifetime
        append!(Ngroup, floor(Int, 330 * (amax * 100)))
        append!(decay, lifetime)
        append!(pk, 0 + 1/3 * lifetime)
        l = sample_lat(diagram, 0.0)
        append!(lat, l)
        append!(phase, sample_phase(l))
    end
    append!(N, sum(Ngroup))

    ti = 0.0
    temp = deepcopy(Ngroup)
    while ti <= t[length(t)]
        ##lengthi = length(N)
        for i in 1:1:length(Ngroup)
            if decay[i] <= (ti + step)
                filter!(x -> x!=Ngroup[i], temp)
        #lengthi = lengthi - 1
            end
        end
        #Ngroup = temp
        #println(Ngroup, temp)

        fai0 = 0.5 - diagram.t0 / diagram.cycle
        appearance_rate = 10 * (-0.5 * cos(2 * pi * (ti / diagram.cycle + fai0)) + 0.5) + 0.5
        p = Poisson(appearance_rate * step)
        deltaN = rand(p)
        #judge how many spot group exist
        if deltaN != 0
            for i in 1:1:deltaN
                lifetime = sample_lifetime(decay_scale)
                amax = 1e-4 * lifetime
                append!(Ngroup, floor(Int, 330 * (amax * 100)))
                append!(temp, floor(Int, 330 * (amax * 100)))
                append!(decay, lifetime)
                append!(pk, ti + 1/3 * lifetime)
                l = sample_lat(diagram, ti)
                append!(lat, l)
                append!(phase, sample_phase(l))
            end
        else
            append!(temp, 0)
        end
        ti += step
    append!(N, sum(temp))
    end
    return N, decay, pk, lat, phase
end

function conf_spotmodel(filename::String, t::Array{Float64,1}; ap_evo=true)
    #=
    (multi spots scenario)
    if considering the below conditions:(Dumusque et al. 2011c)
    1)form in group (same property, geometrical and lifetime)
    2)the appearance of spot group is relavant with the magnetic cycle
    3)the filling factors of spots are strictly linear with their lifetime
    4)spots linearly grow and die down
    5)differential rotation
    =#
    conf = readin(filename)
    nspot = retrieve(conf, "spots", "nspot", Int)
    rstar = retrieve(conf, "star", "radius", Float64)

    vconv = retrieve(conf, "star", "vconv", Float64)

    t0 = rand(t)
    cycle = retrieve(conf, "diagram", "cycle", Float64)
    init_lat = retrieve(conf, "diagram", "init_lat", Float64)
    step = retrieve(conf, "diagram", "step", Float64)
    decay_scale = retrieve(conf, "spots", "decay_scale", Float64)
    diagram = Diagram(t0, cycle, init_lat)

    N, decay, pk, lat, phase = initial_condition(nspot, step, decay_scale, t, diagram)
    amax = decay .* 1e-4
    lat = lat .* DEG2RAD

    Ngroup = length(decay)[1]

    period = ones(Number, Ngroup) .* retrieve(conf, "star", "period", Float64)
    u = ones(Number, Ngroup) .* retrieve(conf, "star", "u", Float64)
    cspot = ones(Number, Ngroup) .* retrieve(conf, "spots", "cspot", Float64)
    cfac = ones(Number, Ngroup) .* retrieve(conf, "spots", "cfac", Float64)
    incl = ones(Number, Ngroup) .* (retrieve(conf, "star", "incl", Float64) * DEG2RAD)
    Q = ones(Number, Ngroup) .* retrieve(conf, "spots", "Q", Float64)

    SpotModel = RotationSpot(nspot, rstar, period, u, cspot, vconv, cfac,
        amax, pk, decay, phase, incl, lat, Q)

    diffrot = retrieve(conf, "spots", "differential", Float64)
    omega_eq = 2 * pi / retrieve(conf, "star", "period", Float64)
    differential = Differential(omega_eq, diffrot)

    return [SpotModel, differential]
end
