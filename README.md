# Diabatizem

```Julia
using Diabatizem
config = loadConfiguration("configuration.json");
data = buildData(config.input_data.hamiltonian_adiabatic, config.input_data.coupling_∂_∂R_adiabatic, config.settings.interpolation);
areas = detectSinglePeakAreas(data.∂_∂R, config.settings.nonadiabatic_areas, 40.0);
lz = detectLandauZenerAreas(data.Hₐ, areas, config.settings.nonadiabatic_areas, 40.0);
lzz = fitSinglePeakCouplings(lz);
pwf = deriveLandauZenerFunctions(lzz);
```

