using Plots

x_ticks=Vector{Int}(); append!(x_ticks, 1:9); append!(x_ticks, 10:5:100)
x_tick_labels=(x_ticks, collect("$tick" for tick in x_ticks))
pl = plot(title="Diabatization Transforming Matrix (NaH Quasimolecule)",
  xscale = :log10,
  xlims=(1, 100), xlabel = "R, Bohr", xticks=x_tick_labels,
  ylims=(-1, 1), ylabel = "S, a.u.", legendfont=font(12, "DejaVu Sans Mono"), size=(1440, 900));
Nˡ = size(Sᵈᵃᵗᵃ, 2)
N = convert(Int, √(Nˡ))
for l = 1:Nˡ
  i, j = mpos(l, N)
  println("plotting $i, $j")
  plot!(pl, Rᵖᵒⁱⁿᵗˢ, Sᵈᵃᵗᵃ[:, l], label="S$(int2indexsub(i))$(int2indexsub(j))");
end
pl
