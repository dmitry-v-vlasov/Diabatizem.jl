using Calculus
using Logging
using ProgressMeter
using DataFrames

import Dierckx

import CubicEquation

#const tanˡⁱᵐ = 57.28996163075943 # tan(89°)
const tanˡⁱᵐ = 5.671281819617709 # tan(80°)
const degree = π / 180.0
const inv_ϕ₀ = (1 - 1 / golden)
const ϕ₀ = 1/golden

function clearGrid(R::Vector{Float64}, ϵᴿ)
    Rᵛ = Vector{Float64}()
    Rˡ = R[1]
    push!(Rᵛ, Rˡ)
    for Rⁱ ∈ R[2:end]
        if abs(Rⁱ - Rˡ) > ϵᴿ
            Rˡ = Rⁱ
            push!(Rᵛ, Rⁱ)
        else
            # skip
        end
    end
    return Rᵛ
end

function sigmoid_of_name(f1::Function, f2::Function, x₀, α, name::AbstractString)
  return x -> begin
    Logging.configure(level=INFO)
    sf = (1 - sigmoid(x, x₀, α))*f1(x) + sigmoid(x, x₀, α)*f2(x)
    if sf === NaN || sf === NaN64 || sf === NaN32 || sf === NaN16
      error("SIGMOID NaN: x=$x, x₀=$x₀, α=$α, f1(x)=$(f1(x)), f2(x)=$(f2(x)), σ(x)=$(sigmoid(x, x₀, α)), name=$name")
    end
    return sf
  end
end
const sigmoid = (x, x₀, α) -> 1/(1+exp(-(x - x₀)/α))

function mat2string(M::Array{Float64, 2})
    io = IOBuffer()
    Base.showarray(io, M, false)
    return takebuf_string(io)
end

function mat2string(M::Array{Function, 2})
    io = IOBuffer()
    Base.showarray(io, M, false)
    return takebuf_string(io)
end

function load_data(file::AbstractString; header=true)
  readtable(
    file,
    header = header, separator = ' ',
    allowcomments = true, commentmark = '#',
    skipblanks = true, encoding = :utf8, normalizenames = false
  )
end

function save_data(data::DataFrame, file::AbstractString; header=true)
  writetable(file, data; separator=' ', quotemark=' ', header=header, nastring="EMPTY")
end

function load_data(file::AbstractString; header=true)
  readtable(
    file,
    header = header, separator = ' ',
    allowcomments = true, commentmark = '#',
    skipblanks = true, encoding = :utf8, normalizenames = false
  )
end

function save_data(data::DataFrame, file::AbstractString; header=true)
  writetable(file, data; separator=' ', quotemark=' ', header=header, nastring="EMPTY")
end

# ---------------------------------
# Miscellaneous Types
# ---------------------------------
immutable IteratorRow{T<:AbstractMatrix}
  A::T
end
Base.start(::IteratorRow)  = 1
Base.next(it::IteratorRow, i) = (it.A[i, :], i+1)
Base.done(it::IteratorRow, i) = i > size(it.A, 1)

# ---------------------------------
# Miscellaneous utility functions
# ---------------------------------
function filterMatrixRows(M::Array{Float64, 2}, predicate::Function)
  Vᵐ = Vector{Vector{Float64}}()
  for row in IteratorRow(M)
    if predicate(row)
      push!(Vᵐ, row)
    end
  end
  L = size(Vᵐ, 1); N = size(Vᵐ[1], 1)
  Mᵐ = Array{Float64, 2}(L, N)
  for l = 1:L
    Mᵐ[l, :] = Vᵐ[l]
  end
  return Mᵐ
end

function highDerivative(f::Function, x::Float64, Δx::Float64)
  return abs((f(x + Δx) - f(x - Δx)) / (2 * Δx)) > tanˡⁱᵐ
end

function dirtyDerivative(f::Function, x::Float64, Δx::Float64)
  return (f(x + Δx) - f(x - Δx)) / (2 * Δx)
end

function limderivative(f::Function, x::Float64, Δx::Float64)
  fd = (f(x + Δx) - f(x - Δx)) / (2 * Δx)
  return abs(fd) <= tanˡⁱᵐ ?  sderivative(f, x) : sign(fd) * tanˡⁱᵐ
end

function sderivative(f::Function, x::Float64)
  return try
    derivative(f, x)
  catch e
    err("Unable to calculate df/dx on function object $f at the point $x")
    throw(e)
  end
end

function ssecond_derivative(f::Function, x::Float64)
  return try
    second_derivative(f, x)
  catch e
    err("Unable to calculate df/dx on function object $f at the point $x")
    throw(e)
  end
end

"""
  Calculates a column number in a data file for a given
  symmetric or antisymmetric (square) matrix with ZEROs in the main diagonal by
  its row an column indices.

  * i - matrix line
  * j - matrix column
  * N - matrix size

  returns: see sescription
  throws:
    * DomainError in case i is equal to j
"""
function dataColumnOfSymetricMatrix(i::Int, j::Int, N::Int)
  if i < j
    return convert(Int, (2*(N - 1) - i) * (i - 1) / 2 + (j - 1))
  elseif i > j
    return convert(Int, (2*(N - 1) - j) * (j - 1) / 2 + (i - 1))
  else
    throw(DomainError("The equal line and column value '$i' is undefined for the (anti)symmetric matrix with zero main diagonal."))
  end
end

"""
  Calculate size of data table made of a symmetric matrix.
"""
function dataSizeOfSymetricMatrix(N::Int)
  return convert(Int, N*(N - 1) / 2)
end

function sizeOfSymmetricUpperMatrix(Nˡ::Int)
    Logging.configure(level=INFO)
    solver = CubicEquation.Solver()
    roots = solver(0.5, -0.5, -Nˡ)
    info("Roots for Nˡ=$Nˡ: $roots")
    @assert all(isreal, roots)
    N = maximum(filter(x -> x > 0, round(Int, real(roots))))
    L = dataSizeOfSymetricMatrix(N)
    @assert(L == Nˡ, "$L≠$Nˡ")
    return N
end

"""
Convert vector element number to a pair of N×N matrix indices.
"""
function mpos(l::Int, N::Int)
  n, r = divrem(l, N)
  return (r == 0) ? n : n + 1, (r == 0) ? N : r
end

"""
Convert a pair of N×N matrix indices to a one-dimentional array index.
"""
function mvec(i::Int, j::Int, N::Int)
  return N*(i - 1) + j
end

function vec2mat!(v::Vector{Float64}, m::Array{Float64, 2})
  N = size(m, 1)
  @assert N == size(m, 2)
  @assert(N*N == size(v, 1), "With N = $N: $(N*N) ≠ $(size(v, 1)).")
  for l = 1:N*N
    i, j = mpos(l, N)
    m[i, j] = v[l]
  end
end

function mat2vec!(m::Array{Float64, 2}, v::Vector{Float64})
  N = size(m, 1)
  @assert N == size(m, 2)
  @assert N*N == size(v, 1)
  for i = 1:N, j = 1:N
    l = mvec(i, j, N)
    v[l] = m[i, j]
  end
end

function matl2matldiag(Mˡ::Vector{Array{Float64, 2}})
  N = size(Mˡ[1], 1)
  L = size(Mˡ, 1)
  Mˡᵈⁱᵃᵍ = Array{Float64, 2}(L, N)
  for i = 1:L
    H_vector = Vector{Float64}(N)
    for k = 1:N
      H_vector[k] = Mˡ[i][k, k]
    end
    Mˡᵈⁱᵃᵍ[i, :] = H_vector
  end
  return Mˡᵈⁱᵃᵍ
end

function matl2matupper(M::Array{Float64, 2})
  @assert size(M, 1) == size(M, 2) "$(size(M, 1))≠$(size(M, 2))"
  N = size(M, 1)
  L = dataSizeOfSymetricMatrix(N)
  ML = Vector{Float64}(L)
  for i=1:N, j=i+1:N
    l = dataColumnOfSymetricMatrix(i, j, N)
    ML[l] = M[i, j]
  end
  return ML
end

function matl2matlupperx(Mˡ::Vector{Array{Float64, 2}})
  N = size(Mˡ[1], 1)
  L = size(Mˡ, 1); Nᴸ = dataSizeOfSymetricMatrix(N)
  Mˡᵈⁱᵃᵍ = Array{Float64, 2}(L, Nᴸ)
  for l = 1:L
    H_vector = Vector{Float64}(Nᴸ)
    for i = 1:N, j = 1:N
      if i < j && i ≠ j
        k = dataColumnOfSymetricMatrix(i, j, N)
        H_vector[k] = Mˡ[l][i, j]
      end
    end
    Mˡᵈⁱᵃᵍ[l, :] = H_vector
  end
  return Mˡᵈⁱᵃᵍ
end

function matlupperx_ddr2matl(M::Array{Float64, 2})
    Logging.configure(level=INFO)
    L = size(M, 1)
    L_N = size(M, 2)
    N = sizeOfSymmetricUpperMatrix(L_N)
    L_N_check = dataSizeOfSymetricMatrix(N)
    @assert L_N_check == L_N
    Mˡ = Vector{Array{Float64, 2}}(L)
    info("N=$N, L=$L, L_N=$L_N")
    for l = 1:L
        M_l = Array{Float64, 2}(N, N)
        for i = 1:N, j = 1:N
            if i == j
                M_l[i, i] = 0.0
            elseif i < j
                k = dataColumnOfSymetricMatrix(i, j, N)
                M_l[i, j] = M[l, k]
                M_l[j, i] = -M[l, k]
            end
        end
        Mˡ[l] = M_l
    end
    return Mˡ
end

function matl2matdata(Mˡ::Vector{Array{Float64, 2}})
  N = size(Mˡ[1], 1)
  L = size(Mˡ, 1); Nᴸ = N*N
  Mˡᵈⁱᵃᵍ = Array{Float64, 2}(L, Nᴸ)
  for l = 1:L
    H_vector = Vector{Float64}(Nᴸ)
    for i = 1:N, j = 1:N
      k = mvec(i, j, N)
      H_vector[k] = Mˡ[l][i, j]
    end
    Mˡᵈⁱᵃᵍ[l, :] = H_vector
  end
  return Mˡᵈⁱᵃᵍ
end

function matdata2matl(data::DataFrame)
  @assert !isempty(data)
  @assert size(data, 2) > 1
  X = collect(data[:,1])
  Nᵖᵒⁱⁿᵗˢ = size(data, 1)
  L = size(data, 2) - 1
  N = round(Int, √L)
  @assert N*N == L
  Y = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  for k = 1:Nᵖᵒⁱⁿᵗˢ
    Yₖ = Array{Float64, 2}(N, N)
    for l = 1:N*N
      i, j = mpos(l, N)
      Yₖ[i, j] = data[k, l + 1]
    end
    Y[k] = Yₖ
  end
  return X, Y
end

function matf2mat(x::Float64, Mᶠ::Array{Function, 2})
  N₁ = size(Mᶠ, 1); N₂ = size(Mᶠ, 2)
  M = Array{Float64, 2}(N₁, N₂)
  for i = 1:N₁, j = 1:N₂
    M[i, j] = Mᶠ[i, j](x)
  end
  return M
end

function matl2matfsl(X::Vector{Float64}, Mˡ::Vector{Array{Float64, 2}}; behaviour::AbstractString="extrapolate")
  L = length(X)
  N = size(Mˡ[1], 1)
  Mᵈᵃᵗᵃ = matl2mdata(Mˡ)
  Mᶠ = Array{Function, 2}(N, N)
  M_spline = Array{Dierckx.Spline1D, 2}(N, N)
  for k = 1:N*N
    Y = Mᵈᵃᵗᵃ[:, k]
    i, j = mpos(k, N)
    @assert(length(X) == length(Y), "length(X) == length(Y), $(length(X)) ≠ $(length(Y))")
    spl = Dierckx.Spline1D(X, Y; w=ones(length(X)), k=1, bc=behaviour, s=0.0)
    Mᶠ[i, j] = R -> Dierckx.evaluate(spl, R)
    M_spline[i, j] = spl
  end
  return Mᶠ, M_spline
end

function matd2vecfsl(X::Vector{Float64}, Mᵈᵃᵗᵃ::Array{Float64, 2})
  L = length(X)
  Nc = size(Mᵈᵃᵗᵃ, 2)
  Mᶠ = Vector{Function}(Nc)
  M_spline = Vector{Dierckx.Spline1D}(Nc)
  for k = 1:Nc
    Y = Mᵈᵃᵗᵃ[:, k]
    spl = Dierckx.Spline1D(X, Y; w=ones(length(X)), k=1, bc="extrapolate", s=0.0)
    Mᶠ[k] = R -> Dierckx.evaluate(spl, R)
    M_spline[k] = spl
  end
  return Mᶠ, M_spline
end

function matDerivative(X::Float64, splines::Array{Dierckx.Spline1D, 2})
  N = size(splines, 1)
  M = Array{Float64, 2}(N, N)
  for i=1:N, j=1:N
    M[i, j] = Dierckx.derivative(splines[i, j], X; nu=1)
  end
  return M
end

function matl2vec(Mˡ::Vector{Array{Float64, 2}}, i, j)
  N = size(Mˡ, 1)
  M = Vector{Float64}(N)
  for l = 1:N
    M[l] = Mˡ[l][i, j]
  end
  return M
end

function matl2mdata(Mˡ::Vector{Array{Float64, 2}})
  L = size(Mˡ, 1)
  N = size(Mˡ[1], 1)
  Mᵈᵃᵗᵃ = Array{Float64, 2}(L, N*N)
  for k = 1:L
    for i = 1:N, j = 1:N
      l = mvec(i, j, N)
      Mᵈᵃᵗᵃ[k, l] = Mˡ[k][i, j]
    end
  end
  return Mᵈᵃᵗᵃ
end

"""
NTRS (NASA Technical Reports Server)
Iott, J.; Haftka, R. T.; Adelman, H. M.
"Selecting step sizes in sensitivity analysis by finite differences"
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850025225.pdf

Not used at the moment....
"""
function Δhₒₗ(ϵₐ, df_dx, d²f_dx²)
  Φ₀ = abs(d²f_dx²) ⋅ (1 + df_dx*df_dx)
  Φ = Φ₀ > 1e-4 ? Φ₀ : 1e-4
  return 2⋅√(ϵₐ / Φ)
end

function Δhₒₚₜ(ΔRₛₜ, df_dx)
  af = abs(df_dx)
  Δh₀ = ΔRₛₜ / (Φᵧ * (af < 1 ? 1 : af))
  return Δh₀ < ΔRₛₜ ? Δh₀ : ΔRₛₜ
end

function Δhₒₚₜʰ(ΔRₛₜ, df_dx, f)
  af = abs(df_dx)*sqrt(1 + abs(f))
  Δh₀ = ΔRₛₜ / (Φᵧ * (af < 1 ? 1 : af))
  return Δh₀ < ΔRₛₜ ? Δh₀ : ΔRₛₜ
end

function splitx(x₀::Float64, x::Float64, f::Function, m₀::Float64, M₀::Float64, Δxₘᵢₙ::Float64, Δxₘₐₓ::Float64)
  n = splitn(x₀, x, f, m₀, M₀, Δxₘᵢₙ, Δxₘₐₓ)
  return collect(x₀:((x - x₀) / n):x)
end

function splitxn(x₀::Float64, x::Float64, n::Int)
  return collect(x₀:((x - x₀) / n):x)
end

function splitn(x₀::Float64, x::Float64, f::Function, m₀::Float64, M₀::Float64, Δxₘᵢₙ::Float64, Δxₘₐₓ::Float64)
  fₓ₀ = f(x₀); fₓ = f(x)
  fₘᵢₙ = min(abs(fₓ₀), abs(fₓ))
  fₘₐₓ = max(abs(fₓ₀), abs(fₓ))

  fʳᵉᶠ = fₘₐₓ
  δf = abs(fₓ - fₓ₀)
  if δf < m₀ || fʳᵉᶠ < m₀
    return abs(x - x₀) > Δxₘₐₓ ?
      ceil(Int, abs(x - x₀)/abs(Δxₘₐₓ)):
      (abs(x - x₀) > Δxₘᵢₙ ?
        ceil(Int, abs(x - x₀)/abs(Δxₘᵢₙ)) :
        1
      )
  end

  δfᵣₑₗ = δf / fʳᵉᶠ
  # if 1 - δfᵣₑₗ > 1e-2
  #   println("High δf: δf=$δf, fₓ₀=$fₓ₀, fₓ=$fₓ; δfᵣₑₗ=$δfᵣₑₗ, n=$(ceil((δf > 1 ? log(δf) : 1)/(1 - δfᵣₑₗ))) at: x₀=$x₀; x=$x")
  # end
  return (1 - δfᵣₑₗ > 1e-2) ?
    ceil(Int, (δf > 1 ? 10.0*log(δf) : 1)/(1 - δfᵣₑₗ)) :
    (abs(x - x₀) > Δxₘₐₓ ?
      ceil(Int, abs(x - x₀)/abs(Δxₘₐₓ)) :
      (abs(x - x₀) > Δxₘᵢₙ ?
        ceil(Int, abs(x - x₀)/abs(Δxₘᵢₙ)) :
        1
      )
    )
end

const INT_SUBSCRIPT = Dict{Int, Char}(
  0 => '₀',
  1 => '₁',
  2 => '₂',
  3 => '₃',
  4 => '₄',
  5 => '₅',
  6 => '₆',
  7 => '₇',
  8 => '₈',
  9 => '₉',
)
function int2indexsub(i::Int)
  @assert i >= 0
  ds = digits(i)[end:-1:1]
  buffer = IOBuffer()
  for d in ds
    print(buffer, INT_SUBSCRIPT[d])
  end
  return takebuf_string(buffer)
end

function int2molstate(i::Int)
  @assert 1 <= i <= length(⚛⚛_STATES) "The condition '1 ≤ $i ≤ $(length(⚛⚛_STATES))' is false."
  return ⚛⚛_STATES[i]
end

function progressCreate(msg::AbstractString, color::Symbol)
  #return Progress(100, 0.05, msg, ProgressMeter.tty_width(bar), color; barglyphs=BarGlyphs(bar))
  return Progress(100, 0.05, msg, ProgressMeter.tty_width(msg), color)
end

function progressCreate(limit::Int, msg::AbstractString, color::Symbol)
  #return Progress(100, 0.05, msg, ProgressMeter.tty_width(bar), color; barglyphs=BarGlyphs(bar))
  return Progress(limit, 0.05, msg, ProgressMeter.tty_width(msg), color)
end

function progress!(p::Progress, counter::Int, showvalues::Any)
  ProgressMeter.update!(p, counter; showvalues = showvalues)
end

function progress!(p::Progress, counter::Int)
  ProgressMeter.update!(p, counter)
end

function finish!(p::Progress)
  ProgressMeter.finish!(p)
end
