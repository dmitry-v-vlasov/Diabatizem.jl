using Plots
x_ticks=Vector{Int}(); append!(x_ticks, 1:9); append!(x_ticks, 10:5:1000)
x_tick_labels=(x_ticks, collect("$tick" for tick in x_ticks))
pl = plot(title="Adiabatic and Diabatic Potentials of NaH Quasimolecule",
  xscale = :log10,
  xlims=(1, 100), xlabel = "R, Bohr", xticks=x_tick_labels,
  ylims=(-0.08, 0.2), ylabel = "Uᴰ, V; a.u.e.",
  size=(1440, 900));
N = size(Hᴬ, 2)

X = 0:0.01:1000;
for i = 1:N
  plot!(pl, X, Hᴬ[i, i], label = "Vᴬ$(int2indexsub(i)) → $(int2molstate(i))");
end
pl
