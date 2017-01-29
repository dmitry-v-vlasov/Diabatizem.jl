using DataFrames
using Interpolations

include("../util/util.jl")

type Data
  hamiltonian::Array{Function, 2}
  ∂_∂R::Array{Function, 2}
  function Data() new() end
end

function buildData(table_Hₐ::DataFrame, table_∂_∂R::DataFrame)
  data = Data()
  buildHₐ!(table_Hₐ, data)
  build_d_dR!(table_∂_∂R, data, numberOfChannels(table_Hₐ))
  return data
end

function buildHₐ!(table_Hₐ::DataFrame, data::Data)
  X = convert(Array{Float64}, table_Hₐ[1])
  N = numberOfChannels(table_Hₐ)
  data.hamiltonian = Array{Function, 2}(N, N)
  for i = 1:N, j = 1:N
    if i == j
      Y = convert(Array{Float64}, table_Hₐ[i + 1])
      itp = interpolate(Y, BSpline(Cubic(Flat())), OnGrid())
      setindex!(data.hamiltonian, R -> itp[R], i, j)
    else
      setindex!(data.hamiltonian, R -> 0, i, j)
    end
  end
end

function build_d_dR!(table_∂_∂R::DataFrame, data::Data, N::Int)
  X = convert(Array{Float64}, table_∂_∂R[1])
  Nc = size(X, 1) - 1
  data.∂_∂R = Array{Function, 2}(N, N)
  for i = 1:N, j = 1:N
    if i ≠ j
      l = dataColumnOfSymetricMatrix(i, j, N)
      Y = convert(Array{Float64}, table_∂_∂R[l + 1])
      itp = interpolate(i > j ? Y : -Y, BSpline(Cubic(Flat())), OnGrid())
      setindex!(data.∂_∂R, R -> itp[R], i, j)
    else
      setindex!(data.∂_∂R, R -> 0, i, j)
    end
  end
end

function numberOfChannels(table_Hₐ::DataFrame)
  size(table_Hₐ, 2) - 1
end
