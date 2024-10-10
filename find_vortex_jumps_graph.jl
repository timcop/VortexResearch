using Revise
using FourierGPE
using VortexDistributions
using JLD2
using Graphs, GraphPlot, GraphMakie

includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("utils.jl")
@load "sol2.jld2" sol;
@load "sim2.jld2" sim;

# t = 34.8
t =35
psi = sol(t);
X = sim.X
x = X[1]; y = X[2]; z = X[3]
@time vorts_x_slices, vorts_y_slices, vorts_z_slices = findVorticesJumps3D(psi, sim);

# vorts_x_zoom, vorts_y_zoom, vorts_z_zoom = zoom_vortex_points(vorts_x_slices, vorts_y_slices, vorts_z_slices, psi, sim, win=1)
vorts_x_zoom, vorts_y_zoom, vorts_z_zoom = zoom_vortex_points_2(vorts_x_slices, vorts_y_slices, vorts_z_slices, psi, sim, win=0, nz=100)


vorts_x = collect(Iterators.flatten(vorts_x_slices))
vorts_y = collect(Iterators.flatten(vorts_y_slices))
vorts_z = collect(Iterators.flatten(vorts_z_slices))

vorts = vcat(vorts_x, vorts_y, vorts_z)
@time adj_mat = vortex_adjacency_matrix(vorts);
g, scene = vortex_graph(vorts, adj_mat, true, false);
display(scene)

vorts_zoom = vcat(vorts_x_zoom, vorts_y_zoom, vorts_z_zoom)
g, scene = vortex_graph(vorts_zoom, adj_mat, true, false);
display(scene)


@time vort_lines, vort_loops, vort_rings = link_graph_vorts(g, repeat_end_of_ring=true)
plot_iso(psi, X, visible=false)
plot_vortex_lines(vorts_zoom, vort_lines, vort_loops, vort_rings, linewidth=10)


