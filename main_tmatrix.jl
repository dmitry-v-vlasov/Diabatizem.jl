using Diabatizem
using Nullables

C = loadConfiguration("configuration.json");
D = buildData(C.input_data.hamiltonian_adiabatic, C.input_data.coupling_∂_∂R_adiabatic, C.settings.interpolation);

Hᴬ = D.Hₐ;
∂_∂Rᴬ = D.∂_∂R;
∂_∂Rᵈᵃᵗᵃ = convert(Matrix{Float64}, C.input_data.coupling_∂_∂R_adiabatic);

Rᵈᵃᵗᵃ = ∂_∂Rᵈᵃᵗᵃ[:, 1];

#Rᵖᵒⁱⁿᵗˢ = transformationMatrix(Hᴬ, ∂_∂Rᴬ, ∂_∂Rᵐᵒᵈᵉˡ, Rᵈᵃᵗᵃ, C.settings.diabatization)
N = size(Hᴬ, 2)
inp = C.input_paths

@assert !isnull(inp.file_transformation_matrix)
@assert isfile(get(inp.file_transformation_matrix))

S_file_data = load_data(get(inp.file_transformation_matrix); header=true)
@info "Loaded file '$(get(inp.file_transformation_matrix))' with data of size $(size(S_file_data))"
Rᵖᵒⁱⁿᵗˢ, S = matdata2matl(S_file_data)

# (Rᵖᵒⁱⁿᵗˢ,
#   Hᴰ, ∂_∂Rᴰ, S) = diabatize(Hᴬ, ∂_∂Rᴬ, Rᵖᵒⁱⁿᵗˢ, false, S, C.settings.diabatization.use_last_transformation_matrix_from);
(Rᵖᵒⁱⁿᵗˢ,
  Hᴰ, ∂_∂Rᴰ, S) = diabatize(Hᴬ, ∂_∂Rᴬ, Rᵖᵒⁱⁿᵗˢ, false, S, Nullable{Float64}());
Uᴰᵈᵃᵗᵃ = matl2matldiag(Hᴰ);
Hᴰᵈᵃᵗᵃ = matl2matlupperx(Hᴰ);
∂_∂Rᴰᵈᵃᵗᵃ = matl2matlupperx(∂_∂Rᴰ);
Sᵈᵃᵗᵃ = matl2matdata(S);

(∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag) = calculate∂²_∂R²(Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ, size(Uᴰᵈᵃᵗᵃ, 2));

saveData(Rᵖᵒⁱⁿᵗˢ, S, Sᵈᵃᵗᵃ, Uᴰᵈᵃᵗᵃ, Hᴰᵈᵃᵗᵃ, ∂_∂Rᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag, C.output_paths)
