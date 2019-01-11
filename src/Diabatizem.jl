module Diabatizem

using JSON
using DataFrames
using DataStructures

include("constants/constants.jl")
include("util/util.jl")
include("util/piecewise.jl")
include("configuration/configuration.jl")
include("calculation/types.jl")
include("data/data.jl")
include("calculation/area.jl")
include("calculation/fit.jl")
include("calculation/equations.jl")
include("calculation/extras.jl")
include("calculation/couplings.jl")
include("calculation/solver.jl")

export Configuration
export NonadiabaticArea, SinglePeakNonadiabaticArea
export LocalSolution

export loadConfiguration
export buildData, saveData

export detectSinglePeakAreas
export detectLandauZenerAreas
export filterSelectedLandauZenerAreas
export fitLandauZenerCouplings
export deriveLandauZenerCouplingFunctions

export solverTransformationMatrixForAreas

export expandLocalSolutions

export loadInitialConditions

export dataColumnOfSymetricMatrix, dataSizeOfSymetricMatrix
export transformationMatrix, error_S
export diabatize
export calculate∂²_∂R²

export matl2matldiag, matl2matlupperx, matl2mdata, matd2vecfsl, matl2matdata, matdata2matl, mat2string
export matf2matl, matl2matupper_data
export mpos, mvec
export int2indexsub, int2molstate
export load_data, save_data
export matlupperx_ddr2matl

export piecewise_function
export findOptimizableCouplings, optimizeCouplings

export saveMatrixElementTable

end # module
