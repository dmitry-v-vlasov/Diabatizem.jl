module Diabatizem

using JSON
using DataFrames
using Logging

export Configuration
export NonadiabaticArea, SinglePeakNonadiabaticArea

export loadConfiguration
export buildData

export detectSinglePeakAreas
export detectLandauZenerAreas
export fitLandauZenerCouplings
export deriveLandauZenerCouplingFunctions

export dataColumnOfSymetricMatrix, dataSizeOfSymetricMatrix
export transformationMatrix, error_S
export diabatize

export matl2matldiag, matl2matlupperx, matl2mdata
export mpos, mvec
export int2indexsub, int2molstate

include("constants/constants.jl")
include("util/util.jl")
include("configuration/configuration.jl")
include("data/data.jl")
include("calculation/types.jl")
include("calculation/area.jl")
include("calculation/fit.jl")
include("calculation/solver.jl")

end # module
