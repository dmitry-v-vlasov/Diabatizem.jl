using Plots
x_ticks=Vector{Int}(undef, 0); append!(x_ticks, 5:7); #append!(x_ticks, 10:5:100)
x_tick_labels=(x_ticks, collect("$tick" for tick in x_ticks))
pl = plot(title="Adiabatic and Diabatic Potentials of CaH Quasimolecule",
  xscale = :log10,
  xlims=(5, 7), xlabel = "R, Bohr", xticks=x_tick_labels,
  #ylims=(0.16, 0.17), ylabel = "Uᴰ, V; a.u.e.",
  ylims=(0.0, 0.1), ylabel = "Uᴰ, V; a.u.e.",
  size=(1440, 900));
N = size(Uᴰᵈᵃᵗᵃ, 2)

for i = 1:N
  plot!(pl, Rᵖᵒⁱⁿᵗˢ, Uᴰᵈᵃᵗᵃ[:, i], label = "Uᴰ$(int2indexsub(i))");
end

X = 5:0.01:7;
for i = 1:N
  plot!(pl, X, Hᴬ[i, i], linestyle=:dot, label = "Vᴬ$(int2indexsub(i)) → $(int2molstate(i))");
end
pl
