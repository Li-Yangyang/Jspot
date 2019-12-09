#include("spot.jl")
#include("config.jl")
#include("evolution.jl")

function simulation(SpotModel::Spot; nper=20, npper=1000, comb=false)
    meanperiod = mean(SpotModel.period)
    SpotModel.pk = rand(Normal(2 * meanperiod, (nper + 4) * meanperiod), SpotModel.nspot)
    tmin = 0.0
    tmax = meanperiod * nper
    t = Array(tmin:meanperiod/npper:tmax)
    N = size(t)[1]
    dF, dRV, BIS = calc(SpotModel,t)
    if comb
        dF = sum(dF;dims=1)
        dRV = sum(dRV;dims=1)
        BIS = sum(BIS;dims=1)
    end
    return dF, dRV, BIS
end

function simulation(SpotModel::Spot, diagram::Diagram; nper=20, npper=1000, comb=false)
    meanperiod = mean(SpotModel.period)
    SpotModel.pk = rand(Normal(2 * meanperiod, (nper + 4) * meanperiod), SpotModel.nspot)
    tmin = 0.0
    tmax = meanperiod * nper
    t = Array(tmin:meanperiod/npper:tmax)
    diagram.t0 = rand(t)
    SpotModel.lat = [corrected_lat(diagram, SpotModel.pk[i]) * DEG2RAD for i in 1:size(SpotModel.pk)[1]]
    N = size(t)[1]
    dF, dRV, BIS = calc(SpotModel,t)
    if comb
        dF = sum(dF;dims=1)
        dRV = sum(dRV;dims=1)
        BIS = sum(BIS;dims=1)
    end
    return dF, dRV, BIS
end

function simulation(SpotModel::Spot, differential::Differential; nper=20, npper=1000, comb=false)
    meanperiod = mean(SpotModel.period)
    SpotModel.pk = rand(Normal(2 * meanperiod, (nper + 4) * meanperiod), SpotModel.nspot)
    tmin = 0.0
    tmax = meanperiod * nper
    t = Array(tmin:meanperiod/npper:tmax)
    SpotModel.lat = [corrected_lat(diagram, SpotModel.pk[i]) * DEG2RAD for i in 1:size(SpotModel.pk)[1]]
    omega = differential.omega_eq .* (1 .- differential.diffrot  .* (sin.(SpotModel.lat)).^2) # quadratic form
    SpotModel.period = 2 .* pi ./ omega
    N = size(t)[1]
    dF, dRV, BIS = calc(SpotModel,t)
    if comb
        dF = sum(dF;dims=1)
        dRV = sum(dRV;dims=1)
        BIS = sum(BIS;dims=1)
    end
    return dF, dRV, BIS
end

function simulation(SpotModel::Spot, diagram::Diagram, differential::Differential; nper=20, npper=1000, comb=false)
    meanperiod = mean(SpotModel.period)
    SpotModel.pk = rand(Normal(2 * meanperiod, (nper + 4) * meanperiod), SpotModel.nspot)
    tmin = 0.0
    tmax = meanperiod * nper
    t = Array(tmin:meanperiod/npper:tmax)
    diagram.t0 = rand(t)
    SpotModel.lat = [corrected_lat(diagram, SpotModel.pk[i]) * DEG2RAD for i in 1:size(SpotModel.pk)[1]]
    omega = differential.omega_eq .* (1 .- differential.diffrot  .* (sin.(SpotModel.lat)).^2) # quadratic form
    SpotModel.period = 2 .* pi ./ omega
    N = size(t)[1]
    dF, dRV, BIS = calc(SpotModel,t)
    if comb
        dF = sum(dF;dims=1)
        dRV = sum(dRV;dims=1)
        BIS = sum(BIS;dims=1)
    end
    return dF, dRV, BIS
end


#function singlespot(filename)
#    spotmodel = conf_spotmodel(filename)
