## Load packages and 
using Revise
using FourierGPE
using VortexDistributions
using JLD2
includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("utils.jl")

## Load simulation data
@load "sol.jld2" sol;
@load "sim.jld2" sim;

t_idx = 33;
X = sim.X;

psi = sol[t_idx];
plot_iso(psi, X, false, true)

## Detection params
N_interp = 8; α = 0.;
N_radius = N_interp;
ϵ1 = epsilonBall2(X, N_radius, 0.0)
# ϵ2 = epsilonBall2(X, N_radius, 0.5)
# ϵ2 = 0.25
ϵ2 = epsilonBall2(X, N_radius-1, 0.0)

## Detection
vorts_3d = find_vortex_points_3d(psi, X, N_interp);
vorts_class = connect_vortex_points_3d(vorts_3d, X, N_interp, ϵ2,  true);
println("Number of vortices: ", length(vorts_class))
## Plotting
plot_iso(psi, X, false, true)
scatterVortsOnIso!(vorts_3d; markersize=ϵ2, transparency=true, alpha=0.05, shininess=3000, color=:blue)

# scatterVortsOnIso!(vorts_3d; markersize=0.05, color=:red)
scatterClassifiedVortices(vorts_class, vorts_3d, X, true, ϵ2)

