using Diabatizem

C = loadConfiguration("configuration.json")
D = buildData(C.input_data.hamiltonian_adiabatic, C.input_data.coupling_∂_∂R_adiabatic, C.settings.interpolation)

Hᴬ = D.Hₐ
∂_∂Rᴬ = D.∂_∂R

A = detectSinglePeakAreas(∂_∂Rᴬ, C.settings.nonadiabatic_areas, 40.0)
Aˡᶻ = detectLandauZenerAreas(Hᴬ, A, C.settings.nonadiabatic_areas, 40.0)
Aˡᶻᶠ = fitLandauZenerCouplings(Aˡᶻ)
∂_∂Rᵐᵒᵈᵉˡ = deriveLandauZenerCouplingFunctions(Aˡᶻᶠ)

(Rᵖᵒⁱⁿᵗˢ,
  S, Sᵈᵃᵗᵃ) = transformationMatrix(Hᴬ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, C.settings.diabatization)
(Rᵖᵒⁱⁿᵗˢ,
  Hᴰ, ∂_∂Rᴰ) = diabatize(Hᴬ, ∂_∂Rᴬ, Rᵖᵒⁱⁿᵗˢ, S)
Uᴰᵈᵃᵗᵃ = matl2matldiag(Hᴰ)
∂_∂Rᴰᵈᵃᵗᵃ = matl2matlupperx(∂_∂Rᴰ)
