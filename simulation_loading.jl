using Revise
using FourierGPE
using VortexDistributions
using JLD2

@load "sol.jld2" sol;
@load "sim.jld2" sim;

# @unpack_Sim sim;

includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("utils.jl")

t_idx = 80;
psi = sol[t_idx];
X = sim.X;

plot_iso(psi, X)

N_interp = 2;

# This finds all vortex points intersecting the planes in 3 directions
@time vorts_3d = find_vortex_points_3d(psi, X, N_interp) 
plot_iso(psi, X)
scatterVortsOnIso!(vorts_3d, 0.05)


# This creates an array of sets of connected vortices unordered
@time vorts_class = connect_vortex_points_3d(vorts_3d, X, 0.1, N_interp, true);
plot_iso(psi, X, false)

scatterClassifiedVortices(vorts_class, vorts_3d, X, false, 0.04)

# This orders the vortices
@time v_sort = sort_classified_vorts_3d(vorts_class, vorts_3d, X); 
plot_iso(psi, X, false)
periodicPlotting(v_sort, X, 10)

# Plotting epsilon balls around each vortex point
x = X[1]; y = X[2]; z = X[3];
dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
α = 0.;
N_interp = 4;
ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)/N_interp

if N_interp <= 2
    ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)
else
    ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)/(N_interp-1)
end
if ϵ < dx/3
    ϵ = dx/3
end

ϵ2 = sqrt(dx^2 + dy^2 + dz^2)/N_interp
vorts_3d = find_vortex_points_3d(psi, X, N_interp);
vorts_class = connect_vortex_points_3d(vorts_3d, X, α, N_interp, true);

ϵ
ϵ2
balls = [GLMakie.GeometryBasics.Sphere(Point3f(v[1:3]), ϵ) for v in vorts_3d]
balls2 = [GLMakie.GeometryBasics.Sphere(Point3f(v[1:3]), ϵ2) for v in vorts_3d]

plot_iso(psi, X, false, false)
scatterVortsOnIso!(vorts_3d, 1)
# for b in balls
#     mesh!(b, color = :blue, transparency=true, alpha=0.4)
# end
# for b in balls2
#     mesh!(b, color = :red, transparency=true, alpha=0.4)
# end

plot_iso(psi, X, false, false) # Plots the 3d axis

# Balls 1
plot_iso(psi, X, false, false)
mesh!(merge(normal_mesh.(balls)), transparency=true, alpha=0.9, color=:white, shininess=2)
scatterVortsOnIso!(vorts_3d, ϵ)

# Balls 2
mesh!(merge(normal_mesh.(balls2)), transparency=true, alpha=0.4)

ϵ