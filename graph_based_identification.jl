# # t = 34.8 3 rings
# using Revise
# using FourierGPE
# using VortexDistributions
# using JLD2
# using BenchmarkTools
# using FLoops
# using Interpolations

# includet("utils_plots.jl") # includet is short hand for include and then track the changes with revise
# includet("utils.jl")
# @load "sol2.jld2" sol;
# @load "sim2.jld2" sim;

# # psi = sol(34.8);
# psi = sol[20];
# plot_iso(psi, sim.X)

# function vortex_points_on_mesh(psi, X, periodic=true)

#     x = X[1]; y = X[2]; z = X[3];
#     x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);
#     dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];

#     vorts_xslice = []
#     vorts_yslice = []
#     vorts_zslice = []

#     vorts = []

#     grid_round = 0.0

#     for xidx in eachindex(x)

#         ψyz_xi = Torus(psi[xidx, :, :], y, z)
#         vyz_xi = findvortices(ψyz_xi; periodic)
#         for v in vyz_xi

#             lower_dist_y = v.xv > 0 ? v.xv % 0.25 : 0.25 - abs(v.xv % 0.25)
#             upper_dist_y = v.xv > 0 ? 0.25 - (v.xv % 0.25) : abs(v.xv % 0.25)
#             lower_dist_z = v.yv > 0 ? v.yv % 0.25 : 0.25 - abs(v.yv % 0.25)
#             upper_dist_z = v.yv > 0 ? 0.25 - (v.yv % 0.25) : abs(v.yv % 0.25)

#             if min(lower_dist_y, upper_dist_y) < grid_round
#                 # collapse point to nearest grid point
#                 if lower_dist_y < upper_dist_y
#                     v.xv = v.xv - lower_dist_y
#                 else
#                     v.xv = v.xv + upper_dist_y
#                 end
#             end

#             if min(lower_dist_z, upper_dist_z) < grid_round
#                 # collapse point to nearest grid point
#                 if lower_dist_z < upper_dist_z
#                     v.yv = v.yv - lower_dist_z
#                 else
#                     v.yv = v.yv + upper_dist_z
#                 end
#             end
#             vx = [x[xidx], v.xv, v.yv]

#             push!(vorts_xslice, vx)
#         end
#     end

#     for yidx in eachindex(y)
#         ψxz_yi = Torus(psi[:, yidx, :], x, z)
#         vxz_yi = findvortices(ψxz_yi; periodic)
#         # push!(vorts_yslice, vxz_yi)
#         for v in vxz_yi

#             lower_dist_x = v.xv > 0 ? v.xv % 0.25 : 0.25 - abs(v.xv % 0.25)
#             upper_dist_x = v.xv > 0 ? 0.25 - (v.xv % 0.25) : abs(v.xv % 0.25)
#             lower_dist_z = v.yv > 0 ? v.yv % 0.25 : 0.25 - abs(v.yv % 0.25)
#             upper_dist_z = v.yv > 0 ? 0.25 - (v.yv % 0.25) : abs(v.yv % 0.25)

#             if min(lower_dist_x, upper_dist_x) < grid_round
#                 # collapse point to nearest grid point
#                 if lower_dist_x < upper_dist_x
#                     v.xv = v.xv - lower_dist_x
#                 else
#                     v.xv = v.xv + upper_dist_x
#                 end
#             end

#             if min(lower_dist_z, upper_dist_z) < grid_round
#                 # collapse point to nearest grid point
#                 if lower_dist_z < upper_dist_z
#                     v.yv = v.yv - lower_dist_z
#                 else
#                     v.yv = v.yv + upper_dist_z
#                 end
#             end

#             vy = [v.xv, y[yidx], v.yv]
#             push!(vorts_yslice, vy)
#         end
#     end

#     for zidx in eachindex(z)
#         ψxy_zi = Torus(psi[:, :, zidx], x, y)
#         vxy_zi =findvortices(ψxy_zi; periodic)
#         # push!(vorts_zslice, vxy_zi)
#         for v in vxy_zi

#             lower_dist_x = v.xv > 0 ? v.xv % 0.25 : 0.25 - abs(v.xv % 0.25)
#             upper_dist_x = v.xv > 0 ? 0.25 - (v.xv % 0.25) : abs(v.xv % 0.25)
#             lower_dist_y = v.yv > 0 ? v.yv % 0.25 : 0.25 - abs(v.yv % 0.25)
#             upper_dist_y = v.yv > 0 ? 0.25 - (v.yv % 0.25) : abs(v.yv % 0.25)

#             if min(lower_dist_x, upper_dist_x) < grid_round
#                 # collapse point to nearest grid point
#                 if lower_dist_x < upper_dist_x
#                     v.xv = v.xv - lower_dist_x
#                 else
#                     v.xv = v.xv + upper_dist_x
#                 end
#             end

#             if min(lower_dist_y, upper_dist_y) < grid_round
#                 # collapse point to nearest grid point
#                 if lower_dist_y < upper_dist_y
#                     v.yv = v.yv - lower_dist_y
#                 else
#                     v.yv = v.yv + upper_dist_y
#                 end
#             end

#             vz = [v.xv, v.yv, z[zidx]]
#             push!(vorts_zslice, vz)
#         end
#     end
#     return vorts_xslice, vorts_yslice, vorts_zslice
# end

# function vortex_points_moving_window(psi, X)
#     x = X[1]; y = X[2]; z = X[3];
#     x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);
#     vorts_xslice = []
#     vorts_yslice = []
#     vorts_zslice = []

#     for i in eachindex(x)
#         for j in 1:length(y)-3
#             for k in 1:length(y)-3
#                 ψ = Torus(psi[i, j:j+3, k:k+3], y[j:j+3], z[k:k+3])
#                 vorts = findvortices(ψ)
#                 for v in vorts
#                     vx = [x[i], v.xv, v.yv]
#                     push!(vorts_xslice, vx)
#                 end
#             end
#         end
#     end

#     for j in eachindex(y)
#         for i in 1:length(x)-3
#             for k in 1:length(z)-3
#                 ψ = Torus(psi[i:i+3, j, k:k+3], x[i:i+3], z[k:k+3])
#                 vorts = findvortices(ψ)
#                 for v in vorts
#                     vy = [v.xv, y[j], v.yv]
#                     push!(vorts_yslice, vy)
#                 end
#             end
#         end
#     end

#     for k in eachindex(z)
#         for i in 1:length(x)-3
#             for j in 1:length(y)-3
#                 ψ = Torus(psi[i:i+3, j:j+3, k], x[i:i+3], y[j:j+3])
#                 vorts = findvortices(ψ)
#                 for v in vorts
#                     vz = [v.xv, v.yv, z[k]]
#                     push!(vorts_zslice, vz)
#                 end
#             end
#         end
#     end
#     return vorts_xslice, vorts_yslice, vorts_zslice
# end

# function vortex_points_on_mesh_interp(psi, X, N = 1)

#     # @assert N <= 16
#     @assert N >= 1


#     x = X[1]; y = X[2]; z = X[3];
#     x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);
#     dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];


#     @assert size(psi)[1] == length(x)
#     @assert size(psi)[2] == length(y)
#     @assert size(psi)[3] == length(z)

#     x_itp = interpolate(X[1], BSpline(Linear()));
#     y_itp = interpolate(X[2], BSpline(Linear()));
#     z_itp = interpolate(X[3], BSpline(Linear()));

#     x_etp = extrapolate(x_itp, Line())
#     y_etp = extrapolate(y_itp, Line())
#     z_etp = extrapolate(z_itp, Line())

#     psi_itp = interpolate(psi, BSpline(Quadratic(Periodic(OnCell()))))
#     # psi_etp = extrapolate(psi_itp, Periodic())

#     x_range = LinRange(1,length(x),N*(length(x)))
#     y_range = LinRange(1,length(y),N*(length(y)))
#     z_range = LinRange(1,length(z),N*(length(z)))

#     x = LinRange(x[1], x[end], length(x));
#     y = LinRange(y[1], y[end], length(y));
#     z = LinRange(z[1], z[end], length(z));

#     ## loop vectorisation, run in parallel 
#     vorts3d = []
#     vorts_xslice = []
#     vorts_yslice = []
#     vorts_zslice = []

#     results_x = [[] for _ in 1:Threads.nthreads()]
#     results_y = [[] for _ in 1:Threads.nthreads()]
#     results_z = [[] for _ in 1:Threads.nthreads()]

#     let z = z, y=y
#         @floop for xidx in x_range
#             vorts_x = vortex_array(findvortices(Torus(psi_itp(xidx, y_range[1]:y_range[end], z_range[1]:z_range[end]), y, z)))
#             for vidx_x in 1:size(vorts_x)[1]
#                 v_x = vorts_x[vidx_x, :]
#                 vx_x = [x_etp(xidx), v_x[1], v_x[2], v_x[3]]
#                 push!(results_x[Threads.threadid()], vx_x)
#             end
#         end
#     end

#     let x=x, z=z
#         @floop for yidx in y_range
#             vorts_y = vortex_array(findvortices(Torus(psi_itp(x_range[1]:x_range[end], yidx, z_range[1]:z_range[end]), x, z)))
#             for vidx_y in 1:size(vorts_y)[1]
#                 v_y = vorts_y[vidx_y, :]
#                 vy_y = [v_y[1], y_etp(yidx), v_y[2], v_y[3]]
#                 push!(results_y[Threads.threadid()], vy_y)
#             end
#         end
#     end

#     let x=x, y=y
#         @floop for zidx in z_range
#             vorts_z = vortex_array(findvortices(Torus(psi_itp(x_range[1]:x_range[end], y_range[1]:y_range[end], zidx), x, y)))
#             for vidx_z in 1:size(vorts_z)[1]
#                 v_z = vorts_z[vidx_z, :]
#                 vz_z = [v_z[1], v_z[2], z_etp(zidx), v_z[3]]
#                 push!(results_z[Threads.threadid()], vz_z)
#             end
#         end
#     end

#     vorts_xslice = reduce(vcat, results_x)
#     vorts_yslice = reduce(vcat, results_y)
#     vorts_zslice = reduce(vcat, results_z)

#     vorts3d = vcat([vorts_xslice, vorts_yslice, vorts_zslice]...);
#     return vorts_xslice, vorts_yslice, vorts_zslice, vorts3d
# end

# vorts_xslice_window, vorts_yslice_window, vorts_zslice_window = vortex_points_moving_window(psi, sim.X)

# vorts_xslice, vorts_yslice, vorts_zslice = vortex_points_on_mesh(psi, sim.X, false)

# vorts_xslice_itp, vorts_yslice_itp, vorts_zslice_itp, vorts3d_itp = vortex_points_on_mesh_interp(psi, sim.X, 1)

# vorts_xslice_window
# vorts_xslice
# vorts_xslice_itp



# # vorts_xslice = vorts_xslice_window
# # vorts_yslice = vorts_yslice_window
# # vorts_zslice = vorts_zslice_window

# vorts_xslice_filt = []; vorts_yslice_filt = []; vorts_zslice_filt = []
# vorts_xslice_temp = copy(vorts_xslice)


# X = sim.X
# x = X[1]; y = X[2]; z = X[3];
# x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);
# dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];

# dx = dx; dy = dy; dz = dz;


# while length(vorts_xslice_temp) > 0
#     vi = pop!(vorts_xslice_temp)

#     vix = vi[1]; viy = vi[2]; viz = vi[3]
#     viy_lower = viy < 0 ? viy - (viy % dy) - dy : viy - (viy % dy);
#     viy_upper = viy < 0 ? viy - (viy % dy) : viy - (viy % dy) + dy;
#     viz_lower = viz < 0 ? viz - (viz % dz) - dz : viz - (viz % dz);
#     viz_upper = viz < 0 ? viz - (viz % dz) : viz - (viz % dz) + dz;

#     vorts_on_face = findall((x -> (x[1] ≈ vix && x[2] >= viy_lower && x[2] <= viy_upper && x[3] >= viz_lower && x[3] <= viz_upper)), vorts_xslice_temp)
#     v_new = mean(vcat(vorts_xslice_temp[vorts_on_face], [vi]))
#     push!(vorts_xslice_filt, v_new)
#     if length(vorts_on_face) > 0
#         vorts_xslice_temp = vorts_xslice_temp[setdiff(1:end, vorts_on_face)]
#     end
# end

# vorts_yslice_temp = copy(vorts_yslice)

# while length(vorts_yslice_temp) > 0
#     vi = pop!(vorts_yslice_temp)

#     vix = vi[1]; viy = vi[2]; viz = vi[3]
#     vix_lower = vix < 0 ? vix - (vix % dx) - dx : vix - (vix % dx);
#     vix_upper = vix < 0 ? vix - (vix % dx) : vix - (vix % dx) + dx;
#     viz_lower = viz < 0 ? viz - (viz % dz) - dz : viz - (viz % dz);
#     viz_upper = viz < 0 ? viz - (viz % dz) : viz - (viz % dz) + dz;

#     vorts_on_face = findall((x -> (x[2] ≈ viy && x[1] >= vix_lower && x[1] <= vix_upper && x[3] >= viz_lower && x[3] <= viz_upper)), vorts_yslice_temp)
#     v_new = mean(vcat(vorts_yslice_temp[vorts_on_face], [vi]))

#     push!(vorts_yslice_filt, v_new)
#     if length(vorts_on_face) > 0
#         vorts_yslice_temp = vorts_yslice_temp[setdiff(1:end, vorts_on_face)]
#     end
# end

# vorts_zslice_temp = copy(vorts_zslice)

# while length(vorts_zslice_temp) > 0
#     vi = pop!(vorts_zslice_temp)

#     vix = vi[1]; viy = vi[2]; viz = vi[3]
#     vix_lower = vix < 0 ? vix - (vix % dx) - dx : vix - (vix % dx);
#     vix_upper = vix < 0 ? vix - (vix % dx) : vix - (vix % dx) + dx;
#     viy_lower = viy < 0 ? viy - (viy % dy) - dy : viy - (viy % dy);
#     viy_upper = viy < 0 ? viy - (viy % dy) : viy - (viy % dy) + dy;

#     vorts_on_face = findall((x -> (x[3] ≈ viz && x[1] >= vix_lower && x[1] <= vix_upper && x[2] >= viy_lower && x[2] <= viy_upper)), vorts_zslice_temp)
#     v_new = mean(vcat(vorts_zslice_temp[vorts_on_face], [vi]))
#     push!(vorts_zslice_filt, v_new)
#     if length(vorts_on_face) > 0
#         vorts_zslice_temp = vorts_zslice_temp[setdiff(1:end, vorts_on_face)]
#     end
# end


# all_vorts = vcat(vorts_xslice, vorts_yslice, vorts_zslice)
# all_vorts_filt = vcat(vorts_xslice_filt, vorts_yslice_filt, vorts_zslice_filt)

# N_vorts = length(all_vorts)
# N_vorts_filt = length(all_vorts_filt)



# adjacency_matrix = zeros(Bool, N_vorts, N_vorts);
# adjacency_matrix_filt = zeros(Bool, N_vorts_filt, N_vorts_filt);

# fudge = 0.

# @time for i1 in eachindex(vorts_xslice)
#     v1 = vorts_xslice[i1]
#     v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
#     v1x_lower = v1x - dx; v1x_upper = v1x + dx
#     v1y_lower = v1y < 0 ? v1y - (v1y % dy) - dy : v1y - (v1y % dy)
#     v1y_upper = v1y < 0 ? v1y - (v1y % dy) : v1y - (v1y % dy) + dy
#     v1z_lower = v1z < 0 ? v1z - (v1z % dz) - dz : v1z - (v1z % dz)
#     v1z_upper = v1z < 0 ? v1z - (v1z % dz) : v1z - (v1z % dz) + dz


#     v1y_lower -= fudge; v1y_upper += fudge
#     v1z_lower -= fudge; v1z_upper += fudge

#     a = findall((x -> (x[2] >= v1y_lower && x[2] <= v1y_upper && x[3] >= v1z_lower && x[3] <= v1z_upper && x[1] >= v1x_lower && x[1] <= v1x_upper)), all_vorts)
#     println("[$i1] -> $a")
#     for i2 in a
#         adjacency_matrix[i1, i2] = true
#     end
# end

# @time for i1 in eachindex(vorts_yslice)
#     v1 = vorts_yslice[i1]
#     v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
#     v1x_lower = v1x < 0 ? v1x - (v1x % dx) - dx : v1x - (v1x % dx)
#     v1x_upper = v1x < 0 ? v1x - (v1x % dx) : v1x - (v1x % dx) + dx
#     v1y_lower = v1y - dy; v1y_upper = v1y + dy
#     v1z_lower = v1z < 0 ? v1z - (v1z % dz) - dz : v1z - (v1z % dz)
#     v1z_upper = v1z < 0 ? v1z - (v1z % dz) : v1z - (v1z % dz) + dz

#     v1x_lower -= fudge; v1x_upper += fudge
#     v1z_lower -= fudge; v1z_upper += fudge

#     a = findall((x -> (x[1] >= v1x_lower && x[1] <= v1x_upper && x[3] >= v1z_lower && x[3] <= v1z_upper && x[2] >= v1y_lower && x[2] <= v1y_upper)), all_vorts)
#     println("[$i1] -> $a")
#     for i2 in a
#         adjacency_matrix[i1 + length(vorts_xslice), i2] = true
#     end
# end

# @time for i1 in eachindex(vorts_zslice)
#     v1 = vorts_zslice[i1]
#     v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
#     v1x_lower = v1x < 0 ? v1x - (v1x % dx) - dx : v1x - (v1x % dx)
#     v1x_upper = v1x < 0 ? v1x - (v1x % dx) : v1x - (v1x % dx) + dx
#     v1y_lower = v1y < 0 ? v1y - (v1y % dy) - dy : v1y - (v1y % dy)
#     v1y_upper = v1y < 0 ? v1y - (v1y % dy) : v1y - (v1y % dy) + dy
#     v1z_lower = v1z - dz; v1z_upper = v1z + dz

#     v1x_lower -= fudge; v1x_upper += fudge
#     v1y_lower -= fudge; v1y_upper += fudge

#     a = findall((x -> (x[1] >= v1x_lower && x[1] <= v1x_upper && x[2] >= v1y_lower && x[2] <= v1y_upper && x[3] >= v1z_lower && x[3] <= v1z_upper)), all_vorts)
#     println("[$i1] -> $a")
#     for i2 in a
#         adjacency_matrix[i1 + length(vorts_xslice) + length(vorts_yslice), i2] = true
#     end
# end

# @time for i1 in eachindex(vorts_xslice_filt)
#     v1 = vorts_xslice_filt[i1]
#     v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
#     v1x_lower = v1x - dx; v1x_upper = v1x + dx
#     v1y_lower = v1y < 0 ? v1y - (v1y % dy) - dy : v1y - (v1y % dy)
#     v1y_upper = v1y < 0 ? v1y - (v1y % dy) : v1y - (v1y % dy) + dy
#     v1z_lower = v1z < 0 ? v1z - (v1z % dz) - dz : v1z - (v1z % dz)
#     v1z_upper = v1z < 0 ? v1z - (v1z % dz) : v1z - (v1z % dz) + dz
#     a = findall((x -> (x[2] >= v1y_lower && x[2] <= v1y_upper && x[3] >= v1z_lower && x[3] <= v1z_upper && x[1] >= v1x_lower && x[1] <= v1x_upper)), all_vorts_filt)
#     println("[$i1] -> $a")
#     for i2 in a
#         adjacency_matrix_filt[i1, i2] = true
#     end
# end


# @time for i1 in eachindex(vorts_yslice_filt)
#     v1 = vorts_yslice_filt[i1]
#     v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
#     v1x_lower = v1x < 0 ? v1x - (v1x % dx) - dx : v1x - (v1x % dx)
#     v1x_upper = v1x < 0 ? v1x - (v1x % dx) : v1x - (v1x % dx) + dx
#     v1y_lower = v1y - dy; v1y_upper = v1y + dy
#     v1z_lower = v1z < 0 ? v1z - (v1z % dz) - dz : v1z - (v1z % dz)
#     v1z_upper = v1z < 0 ? v1z - (v1z % dz) : v1z - (v1z % dz) + dz
#     a = findall((x -> (x[1] >= v1x_lower && x[1] <= v1x_upper && x[3] >= v1z_lower && x[3] <= v1z_upper && x[2] >= v1y_lower && x[2] <= v1y_upper)), all_vorts_filt)
#     println("[$i1] -> $a")
#     for i2 in a
#         adjacency_matrix_filt[i1 + length(vorts_xslice_filt), i2] = true
#     end
# end

# @time for i1 in eachindex(vorts_zslice_filt)
#     v1 = vorts_zslice_filt[i1]
#     v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
#     v1x_lower = v1x < 0 ? v1x - (v1x % dx) - dx : v1x - (v1x % dx)
#     v1x_upper = v1x < 0 ? v1x - (v1x % dx) : v1x - (v1x % dx) + dx
#     v1y_lower = v1y < 0 ? v1y - (v1y % dy) - dy : v1y - (v1y % dy)
#     v1y_upper = v1y < 0 ? v1y - (v1y % dy) : v1y - (v1y % dy) + dy
#     v1z_lower = v1z - dz; v1z_upper = v1z + dz
#     a = findall((x -> (x[1] >= v1x_lower && x[1] <= v1x_upper && x[2] >= v1y_lower && x[2] <= v1y_upper && x[3] >= v1z_lower && x[3] <= v1z_upper)), all_vorts_filt)
#     println("[$i1] -> $a")
#     for i2 in a
#         adjacency_matrix_filt[i1 + length(vorts_xslice_filt) + length(vorts_yslice_filt), i2] = true
#     end
# end



# using LinearAlgebra
# using Graphs
# using SGtSNEpi
# using SNAPDatasets
# using MetaGraphs
# using GraphPlot
# using GraphMakie
# using Makie

# for i in eachindex(all_vorts_filt)
#     adjacency_matrix_filt[i, i] = false
# end

# g = SimpleGraph()
# add_vertices!(g, length(all_vorts_filt))
# for i in 1:length(all_vorts_filt)
#     for j in 1:length(all_vorts_filt)
#         if adjacency_matrix_filt[i, j]
#             add_edge!(g, i, j)
#         end
#     end
# end
# length(all_vorts_filt)
# size(adjacency_matrix_filt)
# all_vorts_filt = [[v[1], v[2], v[3]] for v in all_vorts_filt]
# v_pos = (adj_mat) -> all_vorts_filt

# all_vorts_filt[1]

# plot_iso(psi, sim.X)
# graphplot!(g, layout=v_pos, node_size=10)

# g
# deg_1 = findall(x -> x == 1, degree(g))
# deg_2g = findall(x -> x > 2, degree(g))


# close_vs = [all_vorts_filt[deg_1[7]]]

# # close_vs = vcat(all_vorts_filt[deg_1[5:6]], [all_vorts_filt[deg_1[8]]], [all_vorts_filt[deg_1[10]]])

# meshscatter!([v[1] for v in close_vs], [v[2] for v in close_vs], [v[3] for v in close_vs], color = :red, markersize = 0.01)

# function euclid(v1, v2)
#     return sqrt((v1[1] - v2[1])^2 + (v1[2] - v2[2])^2 + (v1[3] - v2[3])^2)
# end

# min_val = 10
# min_v = []
# for v1 in all_vorts_filt
#     dist = euclid(v1, all_vorts_filt[65])
#     if dist < min_val && dist != 0
#         min_val = dist
#         min_v = v1
#     end
# end


# min_v
# meshscatter!(min_v[1], min_v[2], min_v[3], color = :green, markersize = 0.01)

# ###########

# v1 = min_v
# v1x = v1[1]; v1y = v1[2]; v1z = v1[3]
# v1x_lower = v1x < 0 ? v1x - (v1x % dx) - dx : v1x - (v1x % dx)
# v1x_upper = v1x < 0 ? v1x - (v1x % dx) : v1x - (v1x % dx) + dx
# v1y_lower = v1y < 0 ? v1y - (v1y % dy) - dy : v1y - (v1y % dy)
# v1y_upper = v1y < 0 ? v1y - (v1y % dy) : v1y - (v1y % dy) + dy
# v1z_lower = v1z - dz; v1z_upper = v1z + dz
# a = findall((x -> (x[1] >= v1x_lower && x[1] <= v1x_upper && x[2] >= v1y_lower && x[2] <= v1y_upper && x[3] >= v1z_lower && x[3] <= v1z_upper)), all_vorts_filt)
# deg_1_vort = [all_vorts_filt[deg_1[7]]]
# min_v


# all_vorts_filt[a]

# close_vorts_a = all_vorts_filt[a]
# meshscatter!([v[1] for v in close_vorts_a], [v[2] for v in close_vorts_a], [v[3] for v in close_vorts_a], color = :blue, markersize = 0.01)

# v1
# all_vorts_filt[a]
# all_vorts_filt[deg_2g][1]
# all_vorts_filt[66] == all_vorts_filt[deg_2g][1]
# all_vorts_filt[155] == min_v

# adjacency_matrix_filt[155, 66]

# findall(x -> x == min_v, all_vorts_filt)
# findall(x -> x == 1, adjacency_matrix_filt[:, 66])
# all_vorts_filt[66]
# all_vorts_filt[findall(x -> x == 1, adjacency_matrix_filt[:, 66])]

