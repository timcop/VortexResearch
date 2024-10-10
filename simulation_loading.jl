using Revise
using FourierGPE
using VortexDistributions
using JLD2

@load "sol2.jld2" sol;
@load "sim2.jld2" sim;

# @unpack_Sim sim;

includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("utils.jl")

t_idx = 90;
psi = sol[t_idx];
psi = sol(34.8);
X = sim.X;

plot_iso(psi, X)

N_interp = 1;

# This finds all vortex points intersecting the planes in 3 directions
@time vorts_3d = find_vortex_points_3d(psi, X, N_interp) 
plot_iso(psi, X)
scatterVortsOnIso!(vorts_3d, markersize=0.05)


# This creates an array of sets of connected vortices unordered
@time vorts_class = connect_vortex_points_3d(vorts_3d, X, 0.1, N_interp, true);
plot_iso(psi, X, false)

scatterClassifiedVortices(vorts_class, vorts_3d, X, false, 0.04)

# This orders the vortices
@time v_sort = sort_classified_vorts_3d(vorts_class, vorts_3d, X); 
plot_iso(psi, X, false)
periodicPlotting(v_sort, X, 10)

