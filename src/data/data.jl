using DataFrames
using Interpolations

include("../util/util.jl")
include("data_interpolation.jl")

type Data
  hamiltonian::Array{Function, 2}
  ∂_∂R::Array{Function, 2}
  function Data() end
end

function buildData(table_Hₐ::DataFrame, table_∂_∂R::DataFrame)
  data = Data()
  buildHₐ!(table_Hₐ, data)
  build_d_dR!(table_∂_∂R, data, numberOfChannels(table_Hₐ))
end

function buildHₐ!(table_Hₐ::DataFrame, data::Data)
  X = convert(Array{Float64}, table_Hₐ[1])
  N = numberOfChannels(table_Hₐ)
  data.hamiltonian = Array{Function, 2}(N, N)
  for i = 1:N, j = 1:N
    if i == j
      Y = convert(Array{Float64}, table_Hₐ[i + 1])
      itp = interpolate(Y, BSpline(Cubic(Flat())), OnGrid())
      data.hamiltonian[i][j] = R -> itp[R]
    else
      data.hamiltonian[i][j] = R -> 0
    end
  end
end

function build_d_dR!(table_∂_∂R::DataFrame, N::Int, data::Data)
  X = convert(Array{Float64}, table_∂_∂R[1])
  Nc = size(X) - 1
  data.hamiltonian = Array{Function, 2}(N, N)
  for i = 1:N, j = 1:N
    if i ≠ j
      l = i > j ? (2(N - 1) - i) * (i - 1) + (j - 1) : (2(N - 1) - j) * (j - 1) + (i - 1)
      Y = convert(Array{Float64}, table_∂_∂R[l + 1])
      itp = interpolate(i > j ? Y : -Y, BSpline(Cubic(Flat())), OnGrid())
      data.∂_∂R[i][j] = R -> itp[R]
    else
      data.∂_∂R[i][j] = R -> 0
    end
  end
end

function numberOfChannels(table_Hₐ::DataFrame)
  size(table_Hₐ[1]) - 1
end
