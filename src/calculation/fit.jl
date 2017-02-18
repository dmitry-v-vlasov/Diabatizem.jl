
"""
Make a simple fitting of a given single peak area with
the Lorenz curve which is defined in the Landau-Zener model as
<1|∂/∂R|2> = τ₀ / ((R - R₀)² + 4τ₀²)
"""
function fitSinglePeakCouplings(areas::Array{Array{SinglePeakNonadiabaticArea, 1}, 2})
  N = size(areas, 1)
  M_Α_lz = Array{Array{LandauZenerArea, 1}, 2}(N, N)
  fill!(M_Α_lz, Array{LandauZenerArea, 1}())
  for i = 1:N, j = 1:N
    M_Α_lz[i, j] = Array{LandauZenerArea, 1}()
    ∂_∂R_areas = areas[i, j]
    for Α in ∂_∂R_areas
      Rₐ = Α.coordinate_from; Rᵦ = Α.coordinate_to
      R₀ = Α.coordinate_∂_∂R
      τₕ₀ = Α.value_∂_∂R; τ₀ = 1/(4τₕ₀)

      Α_lz = LandauZenerArea(Α.states, R₀, τ₀, Rₐ, Rᵦ)
      push!(M_Α_lz[i, j], Α_lz)
    end
  end
  return M_Α_lz
end
