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

# Mutate the solution to xspace before saving to jld2
for i in eachindex(sol)
    sol[i] = xspace(sol[i], sim)
end

@save "sol2.jld2" sol
@save "sim2.jld2" sim

plot_iso(sol(34.8), sim.X)

