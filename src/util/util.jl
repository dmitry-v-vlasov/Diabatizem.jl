using Interpolations

# ---------------------------------
# Miscellaneous utility functions
# ---------------------------------
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
  @assert N*N == size(v, 1)
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

function matf2mat(x::Float64, Mᶠ::Array{Function, 2})
  N₁ = size(Mᶠ, 1); N₂ = size(Mᶠ, 2)
  M = Array{Float64, 2}(N₁, N₂)
  for i = 1:N₁, j = 1:N₂
    M[i, j] = Mᶠ[i, j](x)
  end
  return M
end

function matl2matfsl(X::Vector{Float64}, Mˡ::Vector{Array{Float64, 2}})
  L = size(X, 1)
  N = size(Mˡ[1], 1)
  Mᵈᵃᵗᵃ = matl2mdata(Mˡ)
  Mᶠ = Array{Function, 2}(N, N)
  for k = 1:N*N
    Y = Mᵈᵃᵗᵃ[:, k]
    i, j = mpos(k, N)
    itp = interpolate((X,), Y, Gridded(Linear()))
    Mᶠ[i, j] = R -> itp[R]
  end
  return Mᶠ
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
