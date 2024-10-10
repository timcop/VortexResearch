
## Cell1
using GLMakie
using Random
using FourierGPE
using VortexDistributions
using JLD2

includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("utils.jl")

## Load simulation data


sim_vars = SimulationVariables("sol2.jld2", "sim2.jld2");
vortex_plot = VortexPlot(sim_vars);
vortex_plot_observables = VortexPlotObservables(sim_vars, vortex_plot);

PlotVortexPlotObservables(sim_vars, vortex_plot, vortex_plot_observables);
vortex_plot.fig
