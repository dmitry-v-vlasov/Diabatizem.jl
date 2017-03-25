using Plots

∂_∂RC = C.input_data.coupling_∂_∂R_adiabatic
Hᴬ = D.Hₐ
∂_∂Rᴬ = D.∂_∂R

X = 29.633:0.00000001:29.635;

y_ticks=Vector{Int}(); append!(y_ticks, 0:100:400); append!(y_ticks, 500:100:12500)
x_ticks=Vector{Float64}(); append!(x_ticks, 29.633:0.0001:29.635)
pl = plot(title="\"Accurate\" (Q. Chem.) CaH ⟨8|∂/∂R|9⟩ Diabatization",
  xscale = :identity, xlims=(29.633, 29.635), xticks=x_ticks, xlabel = "R, Bohr",
  ylabel = "⟨B|∂/∂R|C⟩, a.u.", yticks=y_ticks,
  size=(1440, 900));
plot!(pl, Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ[:, 50], label = "⟨8|∂/∂R|9⟩, diabatized`");
plot!(pl, X, ∂_∂Rᵐᵒᵈᵉˡ[8, 9], label = "⟨8|∂/∂R|9⟩, Landau-Zener model")
plot!(pl, X, ∂_∂Rᴬ[8, 9], label = "⟨8|∂/∂R|9⟩, original")
