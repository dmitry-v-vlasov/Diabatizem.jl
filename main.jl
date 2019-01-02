using Diabatizem
using Nullables

C = loadConfiguration("configuration.json");
D = buildData(C.input_data.hamiltonian_adiabatic, C.input_data.coupling_∂_∂R_adiabatic, C.settings.interpolation);

Hᴬ = D.Hₐ;
∂_∂Rᴬ = D.∂_∂R;
∂_∂Rᵈᵃᵗᵃ = convert(Matrix{Float64}, C.input_data.coupling_∂_∂R_adiabatic);

A = detectSinglePeakAreas(∂_∂Rᴬ, ∂_∂Rᵈᵃᵗᵃ, Hᴬ, C.settings.nonadiabatic_areas, 1000.0);
Aˡᶻ = detectLandauZenerAreas(Hᴬ, A, C.settings.nonadiabatic_areas, 1000.0);
Aˡᶻᶠ = fitLandauZenerCouplings(Aˡᶻ);
∂_∂Rᵐᵒᵈᵉˡ = deriveLandauZenerCouplingFunctions(Aˡᶻᶠ);
Rᵈᵃᵗᵃ = ∂_∂Rᵈᵃᵗᵃ[:, 1];

#Rᵖᵒⁱⁿᵗˢ = transformationMatrix(Hᴬ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, Rᵈᵃᵗᵃ, C.settings.diabatization)
N = size(Hᴬ, 2)
inp = C.input_paths
S₀ = isnull(inp.file_transformation_matrix_initial) ? Nullable{Matrix{Float64}}() : loadInitialConditions(get(inp.file_transformation_matrix_initial), N)
(Rᵖᵒⁱⁿᵗˢ,
  S, Sᵈᵃᵗᵃ) = transformationMatrix(Hᴬ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, Rᵈᵃᵗᵃ, S₀, C.settings.diabatization);
(ϵ_S, ϵ_Sᵈᵃᵗᵃ) = error_S(S);
(Rᵖᵒⁱⁿᵗˢ,
  Hᴰ, ∂_∂Rᴰ, S) = diabatize(Hᴬ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, Rᵖᵒⁱⁿᵗˢ, true, S, C.settings.diabatization.use_last_transformation_matrix_from);
Uᴰᵈᵃᵗᵃ = matl2matldiag(Hᴰ);
Hᴰᵈᵃᵗᵃ = matl2matlupperx(Hᴰ);
∂_∂Rᴰᵈᵃᵗᵃ = matl2matlupperx(∂_∂Rᴰ);
Sᵈᵃᵗᵃ = matl2matdata(S);

(∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag) = calculate∂²_∂R²(Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ, size(Uᴰᵈᵃᵗᵃ, 2));

saveData(Rᵖᵒⁱⁿᵗˢ, S, Sᵈᵃᵗᵃ, Uᴰᵈᵃᵗᵃ, Hᴰᵈᵃᵗᵃ, ∂_∂Rᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag, C.output_paths)
