using DataFrames

import Dierckx

type Data
  Hₐ::Array{Function, 2}
  ∂_∂R::Array{Function, 2}
  itp_Hₐ::Array{Dierckx.Spline1D, 2}
  itp_∂_∂R::Array{Dierckx.Spline1D, 2}
  function Data() new() end
end

function buildData(table_Hₐ::DataFrame, table_∂_∂R::DataFrame, interpolationSettings::InterpolationSettings)
  data = Data()
  buildHₐ!(table_Hₐ, data, interpolationSettings.hamiltonian)
  build_d_dR!(table_∂_∂R, data, numberOfChannels(table_Hₐ), interpolationSettings.coupling_∂_∂R)
  return data
end

function buildHₐ!(table_Hₐ::DataFrame, data::Data, interpolationType::InterpolationType)
  X = convert(Vector{Float64}, table_Hₐ[1])
  ΔR = X[2] - X[1]
  N = numberOfChannels(table_Hₐ)
  data.Hₐ = Array{Function, 2}(N, N)
  data.itp_Hₐ = Array{Dierckx.Spline1D, 2}(N, N)
  for i = 1:N, j = 1:N
    if i == j
      Y = convert(Vector{Float64}, table_Hₐ[i + 1])
      spl = Dierckx.Spline1D(X, Y; w=ones(length(X)), k=splineDegree(interpolationType), bc="extrapolate", s=0.0)
      setindex!(data.Hₐ, R -> Dierckx.evaluate(spl, R), i, j)
      setindex!(data.itp_Hₐ, spl, i, j)
    else
      setindex!(data.Hₐ, R -> 0, i, j)
    end
  end
end

function build_d_dR!(table_∂_∂R::DataFrame, data::Data, N::Int, interpolationType::InterpolationType)
  X = convert(Vector{Float64}, table_∂_∂R[1])
  ΔR = X[2] - X[1]
  Nc = size(X, 1) - 1
  data.∂_∂R = Array{Function, 2}(N, N)
  data.itp_∂_∂R = Array{Dierckx.Spline1D, 2}(N, N)
  for i = 1:N, j = 1:N
    if i ≠ j
      l = dataColumnOfSymetricMatrix(i, j, N)
      Y = convert(Vector{Float64}, table_∂_∂R[l + 1])

      spl = Dierckx.Spline1D(X, Y; w=ones(length(X)), k=splineDegree(interpolationType), bc="extrapolate", s=0.0)
      setindex!(data.∂_∂R, R -> Dierckx.evaluate(spl, R), i, j)
      setindex!(data.itp_∂_∂R, spl, i, j)
    else
      setindex!(data.∂_∂R, R -> 0, i, j)
    end
  end
end


function splineDegree(interpolationType::InterpolationType)
  if SPLINE_LINEAR::InterpolationType == interpolationType
    return 1
  elseif SPLINE_QUADRATIC::InterpolationType == interpolationType
    return 2
  elseif SPLINE_CUBIC::InterpolationType == interpolationType
    return 3
  end
end

function numberOfChannels(table_Hₐ::DataFrame)
  size(table_Hₐ, 2) - 1
end
