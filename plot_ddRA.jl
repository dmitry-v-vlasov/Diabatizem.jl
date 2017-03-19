using Plots

∂_∂RC = C.input_data.coupling_∂_∂R_adiabatic
Hᴬ = D.Hₐ
∂_∂Rᴬ = D.∂_∂R

X = 29.63:0.00000001:29.64;
pl = plot(title="\"Accurate\" (Q. Chem.) CaH ⟨8|∂/∂R|9⟩",
  xscale = :identity, xlims=(29.63, 29.64), xticks=29.63:0.001:29.64, xlabel = "R, Bohr",
  ylabel = "⟨B|∂/∂R|C⟩, a.u.",
  size=(1440, 900));
plot!(pl, ∂_∂RC[1], ∂_∂RC[51], label = "⟨|∂/∂R|⟩, config");
plot!(pl, X, ∂_∂Rᴬ[8, 9], label = "⟨|∂/∂R|⟩, function")
