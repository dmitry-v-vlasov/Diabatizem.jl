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

## Landau-Zener Coupling Piecewise Function Calculation
```
using Plots
pl = plot(title="Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves", xscale = :identity, xlims=(0, 30), xticks=0:1:30, xlabel = "R, Bohr", ylabel = "⟨B|∂/∂R|C⟩, a.u.", size=(1440, 900));
plot!(pl, X, pwf[3,4], label = "⟨B|∂/∂R|C⟩, Landau-Zener model");
plot!(pl, X, data.∂_∂R[3,4], label = "⟨B|∂/∂R|C⟩, Quantum-chemical data")
```
![Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves](doc/ddrBC_NaH_comparison.png?raw=true "Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves")
