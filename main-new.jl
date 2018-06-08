using Diabatizem

using Logging
Logging.configure(level=INFO)

info("=== Diabatizem ===")
C = loadConfiguration("configuration-new.json");
D = buildData(C.input_data.hamiltonian_adiabatic, C.input_data.coupling_∂_∂R_adiabatic, C.settings.interpolation);

Hᴬ = D.Hₐ;
∂_∂Rᴬ = D.∂_∂R;
∂_∂Rᵈᵃᵗᵃ = convert(Array{Float64, 2}, C.input_data.coupling_∂_∂R_adiabatic);

A = detectSinglePeakAreas(∂_∂Rᴬ, ∂_∂Rᵈᵃᵗᵃ, Hᴬ, C.settings.nonadiabatic_areas, 1000.0);
Aˡᶻ = detectLandauZenerAreas(Hᴬ, A, C.settings.nonadiabatic_areas, 1000.0);

Aˡᶻᶠˡ = filterSelectedLandauZenerAreas(Aˡᶻ, C.settings.diabatization);

# Aˡᶻᶠ = fitLandauZenerCouplings(Aˡᶻᶠˡ);
# ∂_∂Rᵐᵒᵈᵉˡ = deriveLandauZenerCouplingFunctions(Aˡᶻᶠ);
# Rᵈᵃᵗᵃ = ∂_∂Rᵈᵃᵗᵃ[:, 1];
