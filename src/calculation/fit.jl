using Roots
using Logging

import Optim
using LsqFit

"""
Make a simple fitting of a given single peak area with
the Lorenz curve which is defined in the Landau-Zener model as
⟨1|∂/∂R|2⟩ = τ₀ / ((R - R₀)² + 4τ₀²)
Reference:
"СЕЧЕНИЕ ВОЗБУЖДЕНИЯ И МОДЕЛЬ ЛАНДАУ-ЗИНЕРА"
Беляев Андрей Константинович
Известия Российского государственного педагогического университета им. А.И. Герцена
Выпуск № 15 / том 6 / 2006
http://cyberleninka.ru/article/n/sechenie-vozbuzhdeniya-i-model-landau-zinera
"""
function fitLandauZenerCouplings(areas::Array{Vector{SinglePeakNonadiabaticArea}, 2})
    Logging.configure(level=INFO)
  N = size(areas, 1)
  M_Αˡᶻ = Array{Array{LandauZenerArea, 1}, 2}(N, N)
  fill!(M_Αˡᶻ, Array{LandauZenerArea, 1}())
  for i = 1:N, j = 1:N
    M_Αˡᶻ[i, j] = Array{LandauZenerArea, 1}()
    ∂_∂R_areas = areas[i, j]
    for Α in ∂_∂R_areas
      Rₐ = Α.coordinate_from; Rᵦ = Α.coordinate_to
      R₀ = Α.coordinate_∂_∂R
      τₕ₀ = Α.value_∂_∂R; τ₀ = 1/(4τₕ₀)

      info("Curve fitting for\n$Α...")

      model(R::Vector{Float64}, p::Vector{Float64}) = p[2] ./ ((R - p[1]).^2 + 4 * (p[2])^2)
     #  jacobian_model(R::Vector{Float64}, p::Vector{Float64}) = begin
     #    J = Array{Float64}(length(R),length(p))
     #    J[:,1] = 2*p[2]*(R - p[1])./((R - p[1]).^2 + 4 * (p[2])^2).^2    #dmodel/dp[1]
     #    J[:,2] = ((R - p[1]).^2 - 4 * (p[2])^2)./((R - p[1]).^2 + 4 * (p[2])^2).^2  #dmodel/dp[2]
     #    return J
     # end

      Rdata = Α.R_knots
      Fdata = Α.values_∂_∂R
      p₀ = [R₀, τ₀]
      info("Initial parameters: R₀ = $(p₀[1]), τ₀ = $(p₀[2])")
      fit = LsqFit.curve_fit(model, Rdata, Fdata, p₀; maxIter = 10000)
      pₑ = fit.param
      info("Best fit parameters: R₀ = $(pₑ[1]), τ₀ = $(pₑ[2])")
      R₀ = pₑ[1]; τ₀ = pₑ[2]

      # if i == 8 && j == 9
      #     τ₀ = 1/(4*12275)
      #     warn("Hacked $i, $j: τ₀ = $τ₀")
      # end

      Αˡᶻ = LandauZenerArea(Α.states, R₀, τ₀, Rₐ, Rᵦ)
      push!(M_Αˡᶻ[i, j], Αˡᶻ)
    end
  end

  info()
  info("=========== Landau-Zener Couplings Summary ======================")
  for i = 1:N, j=1:N
    if i < j && j - i == 1
      info("*********** Landau-Zener Areas of ⟨$(i)|∂/∂R|$(j)⟩ ***********")
      for Αˡᶻᵛ in M_Αˡᶻ[i, j]
        info(Αˡᶻᵛ)
      end
      info("=======================================================")
    end
  end
  info("==================================================================")

  return M_Αˡᶻ
end

function deriveLandauZenerCouplingFunctions(M_Αˡᶻ::Array{Vector{LandauZenerArea}, 2})
  Logging.configure(level=INFO)

  N = size(M_Αˡᶻ, 1)
  M_∂_∂Rˡᶻ = Array{Function, 2}(N, N)
  fill!(M_∂_∂Rˡᶻ, R -> 0.0)

  info("𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏 Making Landau-Zener Coupling Functions 𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏")
  for i = 1:N, j = 1:N
    areas = M_Αˡᶻ[i, j]
    areas_sorted = sort(areas, alg = InsertionSort, lt = (Α₁, Α₂) -> Α₁.R₀ < Α₂.R₀)
    L = size(areas_sorted, 1)
    functions = Vector{Tuple{Float64, Float64, Function}}()
    breakpoints = Vector{Float64}()
    if L > 1
      info("Piecewise Function for ⟨$(i)|∂/∂R|$(j)⟩; Intervals - $L")
      for k = 1:L-1
        Αₖ = areas_sorted[k]; Αₖ₁ = areas_sorted[k+1];
        Rᵃ = Αₖ.R₀; Rᵇ = Αₖ₁.R₀;
        Δ_∂_∂R(R) = abs(Αₖ.∂_∂R(R)) - abs(Αₖ₁.∂_∂R(R))
        if Δ_∂_∂R(Rᵃ)*Δ_∂_∂R(Rᵇ) > 0
          warn("Cannot establish a breakpoint via root finding in interval [$Rᵃ, $Rᵇ], Δ(∂/∂R(Rᵃ))=$(Δ_∂_∂R(Rᵃ)), Δ(∂_∂R(Rᵇ))=$(Δ_∂_∂R(Rᵇ)); ⟨$(i)|∂/∂R|$(j)⟩$(int2indexsub(k))(Rᵃ)=$(Αₖ.∂_∂R(Rᵃ)), ⟨$(i)|∂/∂R|$(j)⟩$(int2indexsub(k+1))(Rᵇ)=$(Αₖ₁.∂_∂R(Rᵇ))")
          abs_Δ_∂_∂R(R) = abs(Δ_∂_∂R(R))
          result = Optim.optimize(abs_Δ_∂_∂R, Rᵃ, Rᵇ, Optim.Brent())
          Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ = Optim.minimizer(result)
          warn("Found minimum \"distance\" at R=$Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ; ⟨$(i)|∂/∂R|$(j)⟩$(int2indexsub(k))(Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ)=$(Αₖ.∂_∂R(Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ)), ⟨$(i)|∂/∂R|$(j)⟩$(int2indexsub(k+1))(Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ)=$(Αₖ₁.∂_∂R(Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ))")
          push!(breakpoints, Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ)
        else
          Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ = fzero(Δ_∂_∂R, Rᵃ, Rᵇ)
          push!(breakpoints, Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ)
        end
      end
      info("Piecewise Function for ⟨$(i)|∂/∂R|$(j)⟩; Breakpoints: $breakpoints")

      Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ = areas_sorted[1].Rₐ
      for k = 1:L-1
        Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ = breakpoints[k]
        Αₖ = areas_sorted[k]; Αₖ₁ = areas_sorted[k+1];
        push!(functions, (Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ, Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ, Αₖ.∂_∂R))
        Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ = Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ
      end
      push!(functions, (Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ, areas_sorted[L].Rᵦ, areas_sorted[L].∂_∂R))
    elseif L == 1
      info("Single Function for ⟨$(i)|∂/∂R|$(j)⟩")
      push!(functions, (areas_sorted[1].Rₐ, areas_sorted[1].Rᵦ, areas_sorted[1].∂_∂R))
    else
      info("Trivial Zero Function for ⟨$(i)|∂/∂R|$(j)⟩")
      push!(functions, (0.0, Inf, R -> 0.0))
    end

    ∂_∂Rᵖⁱᵉᶜᵉʷⁱˢᵉ(R) = begin
      nᶠ = searchsortedfirst(breakpoints, R)
      ∂_∂Rᵖⁱᵉᶜᵉ = functions[nᶠ][3]
      return ∂_∂Rᵖⁱᵉᶜᵉ(R)
    end

    M_∂_∂Rˡᶻ[i, j] = ∂_∂Rᵖⁱᵉᶜᵉʷⁱˢᵉ
  end
  info("𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏𝔏")

  return M_∂_∂Rˡᶻ
end
