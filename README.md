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
```Julia
using Plots
pl = plot(title="Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves", xscale = :identity, xlims=(0, 30), xticks=0:1:30, xlabel = "R, Bohr", ylabel = "⟨B|∂/∂R|C⟩, a.u.", size=(1440, 900));
plot!(pl, X, pwf[3,4], label = "⟨B|∂/∂R|C⟩, Landau-Zener model");
plot!(pl, X, data.∂_∂R[3,4], label = "⟨B|∂/∂R|C⟩, Quantum-chemical data")
```

* [A. K. Belyaev. Excitation cross sections and the Landau-Zener model. Herzen University Bulletin (Physical Sciences), vol. 6 (15), pp. 213-228 (2006) [in Russian].](http://cyberleninka.ru/article/n/sechenie-vozbuzhdeniya-i-model-landau-zinera)

![Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves](doc/ddrBC_NaH_comparison.png?raw=true "Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves")

## First Test
![](doc/Uᴰ_NaH_V1_V2_V3_V4_pure_Landau_Zener.png?raw=true)

---
[**Atomic Collision Theory Group** ◀ Department of Theoretical Physics and Astronomy of the Herzen State Pedagogical University of Russia, Departement of Optics of the Saint Petersburg State University](http://quasimol.herzen.spb.ru/who-we-are/research-staff)

