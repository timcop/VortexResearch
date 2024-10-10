using Revise
using FourierGPE
using VortexDistributions
using JLD2
using BenchmarkTools
using FLoops
using Interpolations

includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("utils.jl")
@load "sol2.jld2" sol;
@load "sim2.jld2" sim;

psi = sol[60];
plot_iso(psi, sim.X)

X = sim.X
x = X[1]; y = X[2]; z = X[3];

vorts = findvortices(Torus(psi[32, :, :], y, z))


findvortices(Torus(psi[32, 10:11, 41:42], y, z))

(y[10], y[11])
(z[41], z[42])

psi[32, 10:11, 41:42]

y[10:11]
z[41:42]

y_bounds = 9:12
z_bounds = 40:43

y_bounds = y_bounds .+ 1
y_bounds = y_bounds .- 1
z_bounds = z_bounds .+ 1
z_bounds = z_bounds .- 1    

vortex_array(findvortices(Torus(psi[32, y_bounds, z_bounds], y[y_bounds], z[z_bounds])))
findvortices(Torus(psi[32, 9:12, 4:42], y[9:12], z[41:42]))


psi[32, 9:12, 40:43]
vorts_all = []
for i in 3:length(x)-3
    for j in 3:length(x)-3
        for k in 3:length(x)-3
            ψ = Torus(psi[i, j-1:j+2, k-1:k+2], y[j-1:j+2], z[k-1:k+2])
            vorts = vortex_array(findvortices(ψ))
            if length(vorts) > 0
                temp_v = copy(vorts)
                temp_v[:, 1] .= i
                temp_v[:, 2] = vorts[:, 1]
                temp_v[:, 3] = vorts[:, 2]
                push!(vorts_all, (temp_v, (j, k)))
            end
        end 
    end
end

for v in vorts_all
    println(v)
end