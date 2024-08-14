using FourierGPE
using VortexDistributions
using JLD2

@load "sol.jld2" sol
@load "sim.jld2" sim

# @unpack_Sim sim;

include("utils_plots.jl")

t_idx = 10
psi = sol[t_idx]

plot_iso(psi, sim.X)