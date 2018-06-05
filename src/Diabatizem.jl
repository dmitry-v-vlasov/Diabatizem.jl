module Diabatizem

using JSON
using DataFrames
using Logging

export Configuration
export NonadiabaticArea, SinglePeakNonadiabaticArea

export loadConfiguration
export buildData, saveData

export detectSinglePeakAreas
export detectLandauZenerAreas
export fitLandauZenerCouplings
export deriveLandauZenerCouplingFunctions

export loadInitialConditions

export dataColumnOfSymetricMatrix, dataSizeOfSymetricMatrix
export transformationMatrix, error_S
export diabatize
export calculate∂²_∂R²

export matl2matldiag, matl2matlupperx, matl2mdata, matd2vecfsl, matl2matdata, matdata2matl
export mpos, mvec
export int2indexsub, int2molstate

export load_data, save_data

include("constants/constants.jl")
include("util/util.jl")
include("configuration/configuration.jl")
include("data/data.jl")
include("calculation/types.jl")
include("calculation/area.jl")
include("calculation/fit.jl")
include("calculation/solver.jl")
include("calculation/extras.jl")

end # module
