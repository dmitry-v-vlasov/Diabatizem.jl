using Plots

x_ticks=Vector{Int}(); append!(x_ticks, 1:9); append!(x_ticks, 10:5:100)
x_tick_labels=(x_ticks, collect("$tick" for tick in x_ticks))
pl = plot(title="Product of Sᵀ⋅S Transforming Matrix (NaH Quasimolecule)",
  xscale = :log10,
  xlims=(1, 100), xlabel = "R, Bohr", xticks=x_tick_labels,
  ylims=(-0.1, 1.1), ylabel = "Sᵀ⋅S, a.u.", legendfont=font(12, "DejaVu Sans Mono"), size=(1440, 900));
Nˡ = size(Sᵈᵃᵗᵃ, 2)
N = convert(Int, √(Nˡ))
for l = 1:Nˡ
  i, j = mpos(l, N)
  plot!(pl, Rᵖᵒⁱⁿᵗˢ, ϵ_Sᵈᵃᵗᵃ[:, l], label="(Sᵀ⋅S)$(int2indexsub(i))$(int2indexsub(j))");
end
pl
