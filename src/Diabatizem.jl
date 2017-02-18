module Diabatizem

using JSON
using DataFrames

export Configuration
export loadConfiguration
export buildData
export detectSinglePeakAreas
export detectLandauZenerAreas
export fitSinglePeakCouplings

include("constants/constants.jl")
include("util/util.jl")
include("configuration/configuration.jl")
include("data/data.jl")
include("calculation/types.jl")
include("calculation/area.jl")
include("calculation/fit.jl")
include("calculation/solver.jl")

end # module
