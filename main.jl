using Diabatizem
config = loadConfiguration("configuration.json");
data = buildData(config.input_data.hamiltonian_adiabatic, config.input_data.coupling_∂_∂R_adiabatic, config.settings.interpolation);
areas = detectSinglePeakAreas(data.∂_∂R, config.settings.nonadiabatic_areas, 40.0);
lz = detectLandauZenerAreas(data.Hₐ, areas, config.settings.nonadiabatic_areas, 40.0);
lzz = fitSinglePeakCouplings(lz);
pwf = deriveLandauZenerFunctions(lzz);
(Rp, S) = transformationMatrix(data.Hₐ, data.∂_∂R, pwf, config.settings.diabatization);
(Rp, Hd) = diabatize(data.Hₐ, data.∂_∂R, Rp, S);
U = diagonalHᵈvec(Hd);

