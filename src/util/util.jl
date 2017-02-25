# ---------------------------------
# Miscellaneous utility functions
# ---------------------------------
"""
  Calculates a column number in a data file for with a given
  symmetric or antisymmetric (square) matrix with ZEROs in the main diagonal by
  its line an column indices.

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

function matf2mat(x::Float64, Mᶠ::Array{Function, 2})
  N₁ = size(Mᶠ, 1); N₂ = size(Mᶠ, 2)
  M = Array{Float64, 2}(N₁, N₂)
  for i = 1:N₁, j = 1:N₂
    M[i, j] = Mᶠ[i, j](x)
  end
  return M
end

function matl2vec(Mˡ::Array{Array{Function, 2}, 1}, i, j)
  N = size(Mˡ, 1)
  M = Vector{Float64}(N)
  for l = 1:N
    M[l] = Mˡ[l][i, j]
  end
  return M
end

"""
NTRS (NASA Technical Reports Server)
Iott, J.; Haftka, R. T.; Adelman, H. M.
"Selecting step sizes in sensitivity analysis by finite differences"
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850025225.pdf
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
