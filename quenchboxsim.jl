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

# Take a solution at some time t and transform it back to position space
include("utils_plots.jl")

Δt = t[2]-t[1]
# psi = xspace(sol(t[40]), sim);


psi = sol[40]
sim
psi = xspace(psi, sim);

xspace!(psi, sim)

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

