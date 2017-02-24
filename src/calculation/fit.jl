using Roots

"""
Make a simple fitting of a given single peak area with
the Lorenz curve which is defined in the Landau-Zener model as
<1|∂/∂R|2> = τ₀ / ((R - R₀)² + 4τ₀²)
Reference:
"СЕЧЕНИЕ ВОЗБУЖДЕНИЯ И МОДЕЛЬ ЛАНДАУ-ЗИНЕРА"
Беляев Андрей Константинович
Известия Российского государственного педагогического университета им. А.И. Герцена
Выпуск № 15 / том 6 / 2006
http://cyberleninka.ru/article/n/sechenie-vozbuzhdeniya-i-model-landau-zinera
"""
function fitSinglePeakCouplings(areas::Array{Array{SinglePeakNonadiabaticArea, 1}, 2})
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

      Αˡᶻ = LandauZenerArea(Α.states, R₀, τ₀, Rₐ, Rᵦ)
      push!(M_Αˡᶻ[i, j], Αˡᶻ)
    end
  end
  return M_Αˡᶻ
end

function deriveLandauZenerFunctions(M_Αˡᶻ::Array{Array{LandauZenerArea, 1}, 2})
  N = size(M_Αˡᶻ, 1)
  M_∂_∂Rˡᶻ = Array{Function, 2}(N, N)
  fill!(M_∂_∂Rˡᶻ, R -> 0.0)

  for i = 1:N, j = 1:N
    areas = M_Αˡᶻ[i, j]
    areas_sorted = sort(areas, alg = InsertionSort, lt = (Α₁, Α₂) -> Α₁.R₀ < Α₂.R₀)
    L = size(areas_sorted, 1)
    functions = Vector{Tuple{Float64, Float64, Function}}()
    breakpoints = Vector{Float64}()
    if L > 1
      for k = 1:L-1
        Αₖ = areas_sorted[k]; Αₖ₁ = areas_sorted[k+1];
        Rᵃ = Αₖ.R₀; Rᵇ = Αₖ₁.R₀;
        Δ_∂_∂R(R) = abs(Αₖ.∂_∂R(R)) - abs(Αₖ₁.∂_∂R(R))
        Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ = fzero(Δ_∂_∂R, Rᵃ, Rᵇ)
        push!(breakpoints, Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ)
      end

      Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ = areas_sorted[1].Rₐ
      for k = 1:L-1
        Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ = breakpoints[k]
        Αₖ = areas_sorted[k]; Αₖ₁ = areas_sorted[k+1];
        push!(functions, (Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ, Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ, Αₖ.∂_∂R))
        Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ = Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗ
      end
      push!(functions, (Rᵇʳᵉᵃᵏᵖᵒⁱⁿᵗₚᵣₑᵥᵢₒᵤₛ, areas_sorted[L].Rᵦ, areas_sorted[L].∂_∂R))
    elseif L == 1
      push!(functions, (areas_sorted[1].Rₐ, areas_sorted[1].Rᵦ, areas_sorted[1].∂_∂R))
    else
      push!(functions, (0.0, Inf, R -> 0.0))
    end

    ∂_∂Rᵖⁱᵉᶜᵉʷⁱˢᵉ(R) = begin
      nᶠ = searchsortedfirst(breakpoints, R)
      ∂_∂Rᵖⁱᵉᶜᵉ = functions[nᶠ][3]
      return ∂_∂Rᵖⁱᵉᶜᵉ(R)
    end

    M_∂_∂Rˡᶻ[i, j] = ∂_∂Rᵖⁱᵉᶜᵉʷⁱˢᵉ
  end

  return M_∂_∂Rˡᶻ
end
