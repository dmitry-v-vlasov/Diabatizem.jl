using Plots
pl = plot(title="Uᴰ", xscale = :log10, xlims=(1, 100), xlabel = "R, Bohr", ylims=(-0.08, 0.15), ylabel = "Uᴰ, a.u.e", size=(1440, 900));
plot!(pl, Rp, U[1]);
plot!(pl, Rp, U[2]);
plot!(pl, Rp, U[3]);
plot!(pl, Rp, U[4]);

plot!(pl, 0:0.01:100, data.Hₐ[1,1]);
plot!(pl, 0:0.01:100, data.Hₐ[2,2]);

pl
