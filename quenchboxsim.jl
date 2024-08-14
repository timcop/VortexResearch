using FourierGPE
using VortexDistributions
using JLD2

L=(16.,16.,16.);
N=(64,64,64);
sim = Sim(L,N);
@unpack_Sim sim;

## Initialse sim
# parameters
μ = 25.0;
γ = 0.05;
tf = 4/γ;
Nt = 200;
t = LinRange(0.,tf,Nt);

## Run sim
x,y,z = X;
ψi = randn(N)+im*randn(N);
ϕi = kspace(ψi,sim);

@pack_Sim! sim;

## Evolve in k-space
# import FourierGPE
@time sol = runsim(sim); # will take a few minutes to run.

####### Saving the solution #######


@save "sol.jld2" sol

@save "sim.jld2" sim

# Test loading
@load "sol.jld2" sol
@load "sim.jld2" sim

# Take a solution at some time t and transform it back to position space
include("utils_plots.jl")

Δt = t[2]-t[1]
# psi = xspace(sol(t[40]), sim);


psi = sol[40]
sim
psi = xspace(psi, sim);

xspace!(psi, sim)

t[41]
# This will plot the isosurface of the density of the wavefunction at time t=40.

plot_iso(psi, X)

# Parameters for vortex detection

N_interp = 2;

# This finds all vortex points intersecting the planes in 3 directions
@time vorts_3d = find_vortex_points_3d(psi, X, N_interp) 
plot_iso(psi, X, false)
scatterVortsOnIso!(vorts_3d, 0.05)


# This creates an array of sets of connected vortices unordered
@time vorts_class = connect_vortex_points_3d(vorts_3d, X, 0., N_interp, true)
plot_iso(psi, X, false)

scatterClassifiedVortices(vorts_class, vorts_3d, X)

# This orders the vortices
@time v_sort = sort_classified_vorts_3d(vorts_class, vorts_3d, X); 
plot_iso(psi, X, false)
periodicPlotting(v_sort, X, 10)

############ Interactive Plots ############

using GLMakie


psi = xspace(sol(t[end]), sim); # blank plot

scene = plot_iso(psi, X)


fig = scene.figure
ax = scene.axis
# sl = Slider(label = "Time", fig[2, 1], range = t[1]:0.01:t[end], startvalue = t[1])

sl = SliderGrid(
    fig[2, 1],
    (label = "Time", range = t[1]:0.001:t[end], format = "{:.1f}s", startvalue = t[1]))

function density2(psi)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax
    return density
end

psi = lift(sl.sliders[1].value) do t
    xspace(sol(t), sim)
end

d = lift(psi) do psi1
    density2(psi1)
end

cbarPal= :plasma
cmap = get(colorschemes[cbarPal], LinRange(0,1,100));
cmap2 = [(cmap[i], 0.5) for i in 1:100];

volume!(ax, X[1], X[2], X[3], d, algorithm = :iso, isovalue=0.65,isorange=0.075, colormap=cmap2, transparency=true)

vorts_3d = lift(psi) do psi1
    find_vortex_points_3d(psi1, X, N_interp)
end

vorts_filt = lift(vorts_3d) do v
    filter(a->vortInBounds(a, X), v)
end
vorts = lift(vorts_filt) do v
    vorts3DMatrix(v)
end

points = lift(vorts) do v

    Point3f.(v[:, 1], v[:, 2], v[:, 3])

end

function scatterVorts!(points)
    meshscatter!(points, color="blue", markersize=0.05)
end

scatterVorts!(points)
########## Add vortex classification ##########

# vorts_class = connect_vortex_points_3d(vorts_3d, X, 0., N, true)

vorts_class = lift(vorts_3d) do v
    connect_vortex_points_3d(v, X, 0., N_interp, true)
end

