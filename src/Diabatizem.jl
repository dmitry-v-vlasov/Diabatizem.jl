module Diabatizem

using JSON
using DataFrames

export Configuration
export loadConfiguration

include("configuration/configuration.jl")
include("data/data.jl")

end # module
