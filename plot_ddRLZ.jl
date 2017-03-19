using Plots

X = 0:0.001:5;

pl = plot(title="Comparison of a \"Accurate\" (Q. Chem.) ⟨4|∂/∂R|5⟩ Coupling Before and After Diabatization.",
  xscale = :identity, xlims=(0, 5), xticks=0:1:5, xlabel = "R, Bohr",
  ylabel = "⟨4|∂/∂R|5⟩, a.u.",
  size=(1440, 900));
plot!(pl, X, ∂_∂Rᵐᵒᵈᵉˡ[4, 5], label = "⟨4|∂/∂R|5⟩, piecewise");
plot!(pl, X, ∂_∂Rᴬ[4, 5], label = "⟨4|∂/∂R|5⟩, adiabatic (Quantum-chemical data)")
