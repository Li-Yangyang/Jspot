__precompile__()

module Jspot

using ConfParser, Distributions

export conf_spotmodel, simulation

include("spot.jl")
include("evolution.jl")
include("config.jl")
include("simulation.jl")

end # module
