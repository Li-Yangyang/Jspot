include("limbdarkening.jl")
abstract type Spot end

mutable struct RotationSpot<:Spot
    nspot::Int
    rstar::Float64
    period::Array{Number,1}
    u::Array{Number,1}
    cspot::Array{Number,1}
    vconv::Number
    cfac::Array{Number,1}
    amax ::Array{Float64,1}
    pk ::Array{Float64,1}
    decay ::Array{Float64,1}
    phase ::Array{Number,1}
    incl ::Array{Number,1}
    lat ::Array{Number,1}
    Q ::Array{Number,1}
end

function calci(SpotModel::Spot, t::Array{Float64,1}, index::Int)
    amax = SpotModel.amax[index]
    area = ones(size(t)[1]) * amax
    if (SpotModel.pk[index] != 0) * (SpotModel.decay[index] != 0)
        tt = t .- SpotModel.pk[index]
        l = tt .< 0
        area[l] .*= exp.(-tt[l].^2 ./ 10.0 ./ SpotModel.decay[index].^2)
        l = tt .>= 0
        area[l] .*= exp.(-tt[l].^2 ./ 2.0 / SpotModel.decay[index].^2)
    end
    #Fore-shortening
    lon = 2 * pi .* t ./ SpotModel.period[index] .+ SpotModel.phase[index]
    mu = cos(SpotModel.incl[index]) * sin(SpotModel.lat[index]) .+
            sin(SpotModel.incl[index]) * cos(SpotModel.lat[index]) .* cos.(lon)
    #project
    proj = area .* mu
    proj[mu .<0] .= 0
    if SpotModel.u[index] != 0
        N = size(t)[1]
        spot = zeros(N)
        for j in 1:1:N
            spot[j] = dorren_F(SpotModel.u[index], SpotModel.u[index], SpotModel.cspot[index],
            asin(sqrt(area[j])), acos(mu[j]))
        end
    else
        spot = -1.0 .* proj .* SpotModel.cspot[index]
    end
    fac = proj .* SpotModel.Q[index] .* SpotModel.cfac[index] .* (1 .- mu) # This may implement faculae
    dF = copy(spot .+ fac)

    #RV modeling
    veq = 2 * pi * SpotModel.rstar * 6.96e8 / SpotModel.period[index] / 84600.0
    rv_spot = veq * sin(SpotModel.incl[index]) * cos(SpotModel.lat[index]) .* sin.(lon) .* spot
    rv_fac = SpotModel.Q[index] * SpotModel.vconv .* mu .* proj
    dRV = rv_spot .+ rv_fac
    bis = dRV .* cos.(lon)

    return dF, dRV, bis
end

function calc(SpotModel::Spot, t::Array{Float64,1})
    N = size(t)[1]
    M = size(SpotModel.lat)[1]
    dF = Array{Float64}(undef, M, N)
    dRV = Array{Float64}(undef, M, N)
    BIS = Array{Float64}(undef, M, N)
    for i in 1:1:M
        dFi, dRVi, BISi = calci(SpotModel, t, i)
        dF[i,:] = dFi
        dRV[i,:] = dRVi
        BIS[i,:] = BISi
    end
    return dF, dRV, BIS
end
