using Plots

X = 0:0.01:100;

pl = plot(title="Comparison of a \"Accurate\" (Q. Chem.) ⟨4|∂/∂R|5⟩ Coupling Before and After Diabatization.",
  xscale = :identity, xlims=(0, 10), xticks=0:1:10, xlabel = "R, Bohr",
  ylabel = "⟨4|∂/∂R|5⟩, a.u.",
  size=(1440, 900));
plot!(pl, Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ[:, 28], label = "⟨B|∂/∂R|C⟩, diabatized");
plot!(pl, X, ∂_∂Rᴬ[3, 4], label = "⟨B|∂/∂R|C⟩, adiabatic (Quantum-chemical data)")
