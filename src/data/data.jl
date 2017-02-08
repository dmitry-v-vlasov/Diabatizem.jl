using DataFrames
using Interpolations

type Data
  Hₐ::Array{Function, 2}
  ∂_∂R::Array{Function, 2}
  itp_Hₐ::Array{AbstractInterpolation, 2}
  itp_∂_∂R::Array{AbstractInterpolation, 2}
  function Data() new() end
end

function buildData(table_Hₐ::DataFrame, table_∂_∂R::DataFrame, interpolationSettings::InterpolationSettings)
  data = Data()
  buildHₐ!(table_Hₐ, data, interpolationSettings.hamiltonian)
  build_d_dR!(table_∂_∂R, data, numberOfChannels(table_Hₐ), interpolationSettings.coupling_∂_∂R)
  return data
end

function buildHₐ!(table_Hₐ::DataFrame, data::Data, interpolationType::InterpolationType)
  X = convert(Array{Float64}, table_Hₐ[1])
  ΔR = X[2] - X[1]
  N = numberOfChannels(table_Hₐ)
  data.Hₐ = Array{Function, 2}(N, N)
  data.itp_Hₐ = Array{AbstractInterpolation, 2}(N, N)
  for i = 1:N, j = 1:N
    if i == j
      Y = convert(Array{Float64}, table_Hₐ[i + 1])

      interpolationObject = buildInterpolationObject(interpolationType)
      withKnots = isInterpolationWithKnots(interpolationType)

      if withKnots
        itp = interpolate((X,), Y, interpolationObject)
        setindex!(data.Hₐ, R -> itp[R], i, j)
        setindex!(data.itp_Hₐ, itp, i, j)
      else
        itp = interpolate(Y, interpolationObject, OnGrid())
        setindex!(data.Hₐ, R -> itp[R/ΔR + 1], i, j)
        setindex!(data.itp_Hₐ, itp, i, j)
      end
    else
      setindex!(data.Hₐ, R -> 0, i, j)
    end
  end
end

function build_d_dR!(table_∂_∂R::DataFrame, data::Data, N::Int, interpolationType::InterpolationType)
  X = convert(Array{Float64}, table_∂_∂R[1])
  ΔR = X[2] - X[1]
  Nc = size(X, 1) - 1
  data.∂_∂R = Array{Function, 2}(N, N)
  data.itp_∂_∂R = Array{AbstractInterpolation, 2}(N, N)
  for i = 1:N, j = 1:N
    if i ≠ j
      l = dataColumnOfSymetricMatrix(i, j, N)
      Y = convert(Array{Float64}, table_∂_∂R[l + 1])

      interpolationObject = buildInterpolationObject(interpolationType)
      withKnots = isInterpolationWithKnots(interpolationType)

      if withKnots
        itp = interpolate((X,), i < j ? Y : -Y, interpolationObject)
        setindex!(data.∂_∂R, R -> itp[R], i, j)
        setindex!(data.itp_∂_∂R, itp, i, j)
      else
        itp = interpolate(i < j ? Y : -Y, interpolationObject, OnGrid())
        setindex!(data.∂_∂R, R -> itp[R/ΔR + 1], i, j)
        setindex!(data.itp_∂_∂R, itp, i, j)
      end
    else
      setindex!(data.∂_∂R, R -> 0, i, j)
    end
  end
end

function buildInterpolationObject(interpolationType::InterpolationType)
  if SPLINE_LINEAR::InterpolationType == interpolationType
    return Gridded(Linear())
  elseif SPLINE_QUADRATIC::InterpolationType == interpolationType
    return BSpline(Quadratic(Flat()))
  elseif SPLINE_CUBIC::InterpolationType == interpolationType
    return BSpline(Cubic(Flat()))
  end
end

function isInterpolationWithKnots(interpolationType::InterpolationType)
  return SPLINE_LINEAR::InterpolationType == interpolationType
end

function numberOfChannels(table_Hₐ::DataFrame)
  size(table_Hₐ, 2) - 1
end
