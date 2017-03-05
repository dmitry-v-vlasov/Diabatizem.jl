# Diabatizem

```Julia
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
```

## Landau-Zener Coupling Piecewise Function Calculation
```Julia
using Plots

X = 0:0.01:100;

pl = plot(title="Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves",
  xscale = :identity, xlims=(0, 30), xticks=0:2:30, xlabel = "R, Bohr",
  ylabel = "⟨B|∂/∂R|C⟩, a.u.",
  size=(1440, 900));
plot!(pl, X, ∂_∂Rᵐᵒᵈᵉˡ[3,4], label = "⟨B|∂/∂R|C⟩, Landau-Zener model");
plot!(pl, X, ∂_∂Rᴬ[3,4], label = "⟨B|∂/∂R|C⟩, Quantum-chemical data")
```

* [A. K. Belyaev. Excitation cross sections and the Landau-Zener model. Herzen University Bulletin (Physical Sciences), vol. 6 (15), pp. 213-228 (2006) [in Russian].](http://cyberleninka.ru/article/n/sechenie-vozbuzhdeniya-i-model-landau-zinera)

![Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves](doc/ddrBC_NaH_comparison.png?raw=true "Comparison of a More Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling and Landau-Zener Model Curves")

## Non-adiabatic Coupling Diabatization
![Comparison of an Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling Before and After Diabatization](doc/ddrBC_NaH_diabatized.png?raw=true "Comparison of an Accurate (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling Before and After Diabatization")
```Julia
using Plots

X = 0:0.01:100;

pl = plot(title="Comparison of a \"Accurate\" (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling Before and After Diabatization.",
  xscale = :identity, xlims=(0, 30), xticks=0:2:30, xlabel = "R, Bohr",
  ylabel = "⟨B|∂/∂R|C⟩, a.u.",
  size=(1440, 900));
plot!(pl, Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ[:, 6], label = "⟨B|∂/∂R|C⟩, diabatized");
plot!(pl, X, ∂_∂Rᴬ[3, 4], label = "⟨B|∂/∂R|C⟩, adiabatic (Quantum-chemical data)")
```

## Potential Diabatization
```Julia
using Plots
x_ticks=Vector{Int}(); append!(x_ticks, 1:9); append!(x_ticks, 10:5:100)
x_tick_labels=(x_ticks, collect("$tick" for tick in x_ticks))
pl = plot(title="Adiabatic and Diabatic Potentials of NaH Quasimolecule",
  xscale = :log10,
  xlims=(1, 100), xlabel = "R, Bohr", xticks=x_tick_labels,
  ylims=(-0.08, 0.15), ylabel = "Uᴰ, V; a.u.e.",
  size=(1440, 900));
N = size(Uᴰᵈᵃᵗᵃ, 2)

for i = 1:N
  plot!(pl, Rᵖᵒⁱⁿᵗˢ, Uᴰᵈᵃᵗᵃ[:, i], label = "Uᴰ$(int2indexsub(i))");
end

X = 0:0.01:100;
for i = 1:N
  plot!(pl, X, Hᴬ[i, i], label = "Vᴬ$(int2indexsub(i)) → $(int2molstate(i))");
end
pl
```
![](doc/Uᴰ_NaH_V1_V2_V3_V4_pure_Landau_Zener.png?raw=true)

### Transformation Matrix
![](doc/S_transformation_matrix_NaH.png?raw=true)

---
[**Atomic Collision Theory Group** ◀ Department of Theoretical Physics and Astronomy of the Herzen State Pedagogical University of Russia, Department of Optics of the Saint Petersburg State University](http://quasimol.herzen.spb.ru/who-we-are/research-staff)
