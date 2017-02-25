using Calculus
using Sundials

function diabatize(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2}, Rᵖᵒⁱⁿᵗˢ::Array{Float64, 1}, Sˡⁱˢᵗ::Array{Array{Float64, 2}, 1})
  Nᵖᵒⁱⁿᵗˢ = size(Rᵖᵒⁱⁿᵗˢ, 1)
  Hᵈ = Array{Array{Float64, 2}, 1}(Nᵖᵒⁱⁿᵗˢ)
  for i = 1:Nᵖᵒⁱⁿᵗˢ
    R = Rᵖᵒⁱⁿᵗˢ[i]; S = Sˡⁱˢᵗ[i]
    Hᴬ = matf2mat(R, Hₐ)
    Hᴰ = inv(S)*Hᴬ*S
    Hᵈ[i] = Hᴰ
  end
  return Rᵖᵒⁱⁿᵗˢ, Hᵈ
end

function diagonalHᵈ(Hᵈ::Array{Array{Float64, 2}, 1})
  Nᵖᵒⁱⁿᵗˢ = size(Hᵈ, 1)
  Hᵈⁱᵃᵍ = Array{Vector{Float64}, 1}(Nᵖᵒⁱⁿᵗˢ)
  for i = 1:Nᵖᵒⁱⁿᵗˢ
    N = size(Hᵈ[i], 1)
    H_vector = Vector{Float64}(N)
    for k = 1:N
      H_vector[k] = Hᵈ[i][k, k]
    end
    Hᵈⁱᵃᵍ[i] = H_vector
  end
  return Hᵈⁱᵃᵍ
end

function diagonalHᵈvec(Hᵈ::Array{Array{Float64, 2}, 1})
  N = size(Hᵈ[1], 1)
  Nᵖᵒⁱⁿᵗˢ = size(Hᵈ, 1)
  Hᵈⁱᵃᵍ = Vector{Vector{Float64}}(N)
  for i = 1:N
    H_vector = Vector{Float64}(Nᵖᵒⁱⁿᵗˢ)
    for k = 1:Nᵖᵒⁱⁿᵗˢ
      H_vector[k] = Hᵈ[k][i, i]
    end
    Hᵈⁱᵃᵍ[i] = H_vector
  end
  return Hᵈⁱᵃᵍ
end

function transformationMatrix(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2}, config::DiabatizationSettings)
  Rᵇᵉᵍⁱⁿ = config.coordinate_start; Rᵉⁿᵈ = config.coordinate_stop
  ΔRᵐᵃˣ = config.coordinate_step

  N = size(Hₐ, 1); Nᶜ = dataColumnOfSymetricMatrix(N-1, N, N)
  Rᵖᵒⁱⁿᵗˢ = Vector{Float64}(); S₀ = eye(N, N)

  ΔRᵗᵐᵖ = Vector{Float64}(Nᶜ)
  ΔR_∂_∂R(i::Int, j::Int, R::Float64) = Δhₒₚₜʰ(ΔRᵐᵃˣ, derivative(∂_∂R[i, j], R), ∂_∂R[i, j](R))
  ΔRᵐⁱⁿ(R::Float64, ΔRᵗᵉᵐᵖ::Vector{Float64}) = begin
    for i = 1:N, j = 1:N
      if i < j
        l = dataColumnOfSymetricMatrix(i, j, N)
        ΔRᵗᵉᵐᵖ[l] = ΔR_∂_∂R(i, j, R)
      end
    end
    return minimum(ΔRᵗᵉᵐᵖ)
  end

  R = Rᵇᵉᵍⁱⁿ
  while R <= Rᵉⁿᵈ
    push!(Rᵖᵒⁱⁿᵗˢ, R); R += ΔRᵐⁱⁿ(R, ΔRᵗᵐᵖ)
  end

  S = problemCauchy(
    Rᵖᵒⁱⁿᵗˢ, S₀;
    prod_function = diabatizationODE_function,
    data = (N, ∂_∂R, ∂_∂Rᵐᵒᵈᵉˡ),
    ϵʳᵉˡ = 1e-5, ϵᵃᵇˢ = 1e-10
  )

  return Rᵖᵒⁱⁿᵗˢ, S
end

let
  global diabatizationODE_function
  S = Array{Float64, 2}(); dS_dR = Array{Float64, 2}()
  ∂_∂R = Array{Float64, 2}(); ∂_∂Rᵐᵒᵈᵉˡ = Array{Float64, 2}()

  """
  We do not aspite to have an optimal implementation and
  we prefer to have this function with vivid matrix formulae.
  """
  function diabatizationODE_function(R::Float64, S_v::Vector{Float64}, dS_dR_v::Vector{Float64}, data::Tuple{Int, Array{Function, 2}, Array{Function, 2}})
    F_∂_∂R = data[2]; F_∂_∂Rᵐᵒᵈᵉˡ = data[3]

    N = data[1]
    if size(S, 1) < N || size(dS_dR, 1) < N || size(∂_∂R, 1) < N || size(∂_∂Rᵐᵒᵈᵉˡ, 1) < N
      S = Array{Float64, 2}(N, N); dS_dR = Array{Float64, 2}(N, N)
      ∂_∂R = Array{Float64, 2}(N, N); ∂_∂Rᵐᵒᵈᵉˡ = Array{Float64, 2}(N, N)
    end

    vec2mat!(S_v, S)
    for i = 1:N, j = 1:N
      ∂_∂R[i, j] = F_∂_∂R[i, j](R)
      ∂_∂Rᵐᵒᵈᵉˡ[i, j] = F_∂_∂Rᵐᵒᵈᵉˡ[i, j](R)
    end

    dS_dR = (S*∂_∂R - ∂_∂R*S) - S*∂_∂Rᵐᵒᵈᵉˡ
    # dS_dR = -∂_∂Rᵐᵒᵈᵉˡ*S
    mat2vec!(dS_dR, dS_dR_v)

    return Sundials.CV_SUCCESS
  end
end

function problemCauchy(
  Xᵖᵒⁱⁿᵗˢ::Vector{Float64},
  Y₀::Array{Float64, 2};
  prod_function::Function = nothing,
  data::Tuple{Int, Array{Function, 2}, Array{Function, 2}} = nothing,
  ϵʳᵉˡ::Float64 = 1e-3,
  ϵᵃᵇˢ::Float64 = 1e-6)

  N = size(Y₀, 1)

  Yⁱⁿⁱᵗ = Vector{Float64}(N*N); fill!(Yⁱⁿⁱᵗ, 0)
  for i = 1:N, j = 1:N
    l = mvec(i, j, N); Yⁱⁿⁱᵗ[l] = Y₀[i, j]
  end

  Yʳᵉˢ = Sundials.cvode(
    (Xᵖ, Yᵛ, dYᵛ_dx) -> prod_function(Xᵖ, Yᵛ, dYᵛ_dx, data),
    Yⁱⁿⁱᵗ, Xᵖᵒⁱⁿᵗˢ, nothing;
    integrator = :BDF, reltol = ϵʳᵉˡ, abstol = ϵᵃᵇˢ)

  Nᵖᵒⁱⁿᵗˢ = size(Yʳᵉˢ, 1)
  Nˢᵒˡᵘᵗⁱᵒⁿˢ = size(Yʳᵉˢ, 2)
  if Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ N*N throw(ErrorException("Unexpected: Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ N×N: $Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ $N×$N ($(N*N))")) end
  Y = Array{Array{Float64, 2}, 1}(Nᵖᵒⁱⁿᵗˢ)
  for k = 1:Nᵖᵒⁱⁿᵗˢ
    Yₖ = Array{Float64, 2}(N, N)
    for l = 1:N*N
      i, j = mpos(l, N)
      Yₖ[i, j] = Yʳᵉˢ[k, l]
    end
    Y[k] = Yₖ
  end

  return Y
end
