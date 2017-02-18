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
