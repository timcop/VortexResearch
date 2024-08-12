using FourierGPE
using VortexDistributions

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
@time sol = runsim(sim); # will take a few minutes to run.

# To get the density, we need to transform back to position space

## Transform back to position space

# Take a solution at some time t and transform it back to position space

psi = kspace(sol[1], sim);

