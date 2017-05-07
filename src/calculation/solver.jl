using Calculus
using Sundials
using Logging
using ProgressMeter

import Dierckx

function diabatize(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2}, Rᵖᵒⁱⁿᵗˢ::Vector{Float64}, Sˡ::Vector{Array{Float64, 2}})
  Nᵖᵒⁱⁿᵗˢ = size(Rᵖᵒⁱⁿᵗˢ, 1)
  Hᵈ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  ∂_∂Rᵈ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  Sᶠᵘⁿᶜ, S_spline = matl2matfsl(Rᵖᵒⁱⁿᵗˢ, Sˡ)
  for i = 1:Nᵖᵒⁱⁿᵗˢ
    R = Rᵖᵒⁱⁿᵗˢ[i];
    S = Sˡ[i];
    S⁻¹ = S';
    ∇S = Dierckx.derivative.(S_spline, R; nu=1)
    #∇S = matDerivative(R, S_spline)
    #∇S = Calculus.derivative.(Sᶠᵘⁿᶜ, R) # Calculus, what the f***???!!
    #∇S = dirtyDerivative.(Sᶠᵘⁿᶜ, R, 1e-6)
    Hᴬ = matf2mat(R, Hₐ); ∂_∂Rᴬ = matf2mat(R, ∂_∂R); ∂_∂Rᴹ = matf2mat(R, ∂_∂Rᵐᵒᵈᵉˡ)

    Hᴰ = S⁻¹*Hᴬ*S
    ∂_∂Rᴰ = S⁻¹*∂_∂Rᴬ*S⁻¹ + S⁻¹*∇S
    #∂_∂Rᴰ = ∂_∂Rᴬ - ∂_∂Rᴹ
    #∂_∂Rᴰ = ∂_∂Rᴬ

    Hᵈ[i] = Hᴰ; ∂_∂Rᵈ[i] = ∂_∂Rᴰ
  end
  return Rᵖᵒⁱⁿᵗˢ, Hᵈ, ∂_∂Rᵈ
end

function transformationMatrix(Hₐ::Array{Function, 2},
  ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2},
  Rᵛ::Vector{Float64},
  S₀ᵒʷⁿ::Nullable{Array{Float64, 2}},
  C::DiabatizationSettings)

  Logging.configure(level=INFO)

  # -----------
  N = size(Hₐ, 1)
  ℲR⃛ᵒᵖᵗ=(R₀::Float64, R::Float64, m₀, M₀, ΔRₘᵢₙ, ΔRₘₐₓ)->begin
    n = 1
    for i = 1:N, j = 1:N
      if i < j && j - i == 1
        n = max(splitn(R₀, R, ∂_∂R[i, j], m₀, M₀, ΔRₘᵢₙ, ΔRₘₐₓ), n)
      end
    end
    return splitxn(R₀, R, n)
  end
  # -----------

  Rˢᵗᵃʳᵗ = C.coordinate_start
  Rˢᵗᵒᵖ = C.coordinate_stop
  σR = sign(Rˢᵗᵒᵖ - Rˢᵗᵃʳᵗ)
  ΔRₘᵢₙ = σR * min(abs(C.coordinate_step[1]), abs(C.coordinate_step[2]))
  ΔRₘₐₓ = σR * max(abs(C.coordinate_step[1]), abs(C.coordinate_step[2]))

  # -----------
  progress = progressCreate(100, "Making R grid: ", :yellow); progressᶜ = 0
  Rᵖᵒⁱⁿᵗˢ = Vector{Float64}()
  Rₘᵢₙ = min(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ); Rₘₐₓ = max(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ)
  R⃜ = filter(R -> Rₘᵢₙ<=R<=Rₘₐₓ, σR > 0 ? Rᵛ : reverse(unique(Rᵛ))); L = length(R⃜)
  for l = 1:L-1
    Rˡ = R⃜[l]; Rˡ⁺¹ = R⃜[l + 1]
    R⃛ˡ = ℲR⃛ᵒᵖᵗ(Rˡ, Rˡ⁺¹, 2e-2, 1.2e4, ΔRₘᵢₙ, ΔRₘₐₓ)
    for k = 1:length(R⃛ˡ) - 1
      Rˡₖ = R⃛ˡ[k]
      push!(Rᵖᵒⁱⁿᵗˢ, Rˡₖ)
      progressᶜ = floor(Int, (1 - floor(abs(Rˡₖ - Rˢᵗᵒᵖ)/Rₘₐₓ)) * 100)
      progress!(progress, progressᶜ, [(:R, Rˡₖ), (:ΔR, abs(Rˡₖ - R⃛ˡ[k + 1]))])
    end
  end
  push!(Rᵖᵒⁱⁿᵗˢ, R⃜[L])
  finish!(progress)
  # -----------


  S₀ = isnull(S₀ᵒʷⁿ) ? eye(N, N) : get(S₀ᵒʷⁿ)
  S, Sᵈᵃᵗᵃ = problemCauchy(
    Rᵖᵒⁱⁿᵗˢ, S₀;
    prod_function = diabatizationODE_function,
    data = (N, ∂_∂R, ∂_∂Rᵐᵒᵈᵉˡ),
    ϵʳᵉˡ = 1e-5, ϵᵃᵇˢ = 1e-10
  )

  return Rᵖᵒⁱⁿᵗˢ[end:-1:1], S[end:-1:1], Sᵈᵃᵗᵃ[end:-1:1, 1:1:end]
end

function error_S(S::Vector{Array{Float64, 2}})
  L = size(S, 1)
  N = size(S[1], 1)
  ϵ_S = Vector{Array{Float64, 2}}(L)
  for l = 1:L
    ϵ_S[l] = S[l]'*S[l]
  end
  ϵ_Sᵈᵃᵗᵃ = matl2mdata(ϵ_S)
  return ϵ_S, ϵ_Sᵈᵃᵗᵃ
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

    #dS_dR = (S*∂_∂R - ∂_∂R*S) - S*∂_∂Rᵐᵒᵈᵉˡ
    dS_dR = -∂_∂Rᵐᵒᵈᵉˡ*S

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
  Y = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  for k = 1:Nᵖᵒⁱⁿᵗˢ
    Yₖ = Array{Float64, 2}(N, N)
    for l = 1:N*N
      i, j = mpos(l, N)
      Yₖ[i, j] = Yʳᵉˢ[k, l]
    end
    Y[k] = Yₖ
  end

  return Y, Yʳᵉˢ
end
