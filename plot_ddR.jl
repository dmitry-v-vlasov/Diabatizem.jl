using Plots

X = 0:0.01:100;

pl = plot(title="Comparison of a \"Accurate\" (Q. Chem.) NaH ⟨B|∂/∂R|C⟩ Coupling Before and After Diabatization.",
  xscale = :identity, xlims=(0, 30), xticks=0:2:30, xlabel = "R, Bohr",
  ylabel = "⟨B|∂/∂R|C⟩, a.u.",
  size=(1440, 900));
plot!(pl, Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ[:, 6], label = "⟨B|∂/∂R|C⟩, diabatized");
plot!(pl, X, ∂_∂Rᴬ[3, 4], label = "⟨B|∂/∂R|C⟩, adiabatic (Quantum-chemical data)")
