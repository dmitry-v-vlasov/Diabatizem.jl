using Diabatizem

@info "=== Diabatizem ==="
C = loadConfiguration("configuration-new.json");
D = buildData(C.input_data.hamiltonian_adiabatic, C.input_data.coupling_∂_∂R_adiabatic, C.settings.interpolation);

Hᴬ = D.Hₐ;
∂_∂Rᴬ = D.∂_∂R;
∂_∂Rᵈᵃᵗᵃ = convert(Matrix{Float64}, C.input_data.coupling_∂_∂R_adiabatic);

A = detectSinglePeakAreas(∂_∂Rᴬ, ∂_∂Rᵈᵃᵗᵃ, Hᴬ, C.settings.nonadiabatic_areas, 1000.0);
Aˡᶻ = detectLandauZenerAreas(Hᴬ, A, C.settings.nonadiabatic_areas, 1000.0);

Aˡᶻᶠˡ = filterSelectedLandauZenerAreas(Aˡᶻ, C.settings.diabatization);

Aˡᶻᶠ = fitLandauZenerCouplings(Aˡᶻᶠˡ);
∂_∂Rᵐᵒᵈᵉˡ = deriveLandauZenerCouplingFunctions(Aˡᶻᶠ);
Rᵈᵃᵗᵃ = ∂_∂Rᵈᵃᵗᵃ[:, 1];

Sl = solverTransformationMatrixForAreas(Aˡᶻᶠˡ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, Rᵈᵃᵗᵃ, C.settings.diabatization);

expandLocalSolutions(Sl, Rᵈᵃᵗᵃ, C.settings.diabatization);
Rᶜ, S, Hᴰ, ∂_∂Rᴰ, ∂_∂Rᵐ = diabatize(Hᴬ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, Rᵈᵃᵗᵃ, Sl, C.settings.diabatization, Aˡᶻᶠ);

@info "Transformed in the interval [$(Rᶜ[1]), $(Rᶜ[end])]"
@assert !isempty(S)
@info "Transformation matrix size: $(size(S))"
@assert !isempty(Hᴰ)
@info "Hamiltonian matrix size: $(size(Hᴰ))"
@assert !isempty(∂_∂Rᴰ)
@info "⟨·|∂/∂Rᴰ|·⟩ matrix size: $(size(∂_∂Rᴰ))"

@info "type - $(typeof(Hᴰ))"

@info "Matrix data..."
Uᴰᵈᵃᵗᵃ = matl2matldiag(Hᴰ);
Hᴰᵈᵃᵗᵃ = matl2matlupperx(Hᴰ);
∂_∂Rᴰᵈᵃᵗᵃ = matl2matlupperx(∂_∂Rᴰ);
∂_∂Rᵐᵈᵃᵗᵃ = matl2matlupperx(∂_∂Rᵐ);
Sᵈᵃᵗᵃ = matl2matdata(S);
@info "Done."

@info "⟨·|∂²_∂R²ᴰ|·⟩..."
(∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag) = calculate∂²_∂R²(Rᵈᵃᵗᵃ, ∂_∂Rᴰᵈᵃᵗᵃ, size(Uᴰᵈᵃᵗᵃ, 2));
@info "Done."

@info "Saving..."
saveData(Rᶜ, S, Sᵈᵃᵗᵃ, Uᴰᵈᵃᵗᵃ, Hᴰᵈᵃᵗᵃ,
    ∂_∂Rᴰᵈᵃᵗᵃ, ∂_∂Rᵐᵈᵃᵗᵃ, Rᵈᵃᵗᵃ,
    ∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag, Sl, C.output_paths)
@info "Done."
