function epsilonBall1(X, N_interp, α)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)/N_interp

    if N_interp <= 2
        ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)
    else
        ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)/(N_interp-1)
    end
    if ϵ < dx/3
        ϵ = dx/3
    end
    return ϵ
end

function epsilonBall2(X, N_interp, α)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    ϵ = (1+α)*sqrt(dx^2 + dy^2 + dz^2)/N_interp
    return ϵ
end

# Loading sim function

# Interactive plotting functions
function load_sim(sol_path, sim_path)
    sol = JLD2.load(sol_path, "sol")
    sim = JLD2.load(sim_path, "sim")
    X = sim.X
    t = sim.t
    return sol, sim, X, t
end

struct SimulationVariables
    sol::Any
    sim::Any
    X::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    t::Any
end

function SimulationVariables(sol_path, sim_path)
    sol, sim, X, t = load_sim(sol_path, sim_path)
    return SimulationVariables(sol, sim, X, t)
end

mutable struct VortexPlot
    fig::Figure
    toggles::Vector{Toggle}
    toggle_labels::Vector{Label}
    sliders::Vector{Slider}
    interp_depth::Observable
end

function VortexPlot(sim_vars)
    t = sim_vars.t
    fig = Figure()
    lscene = LScene(fig[1, 1])

    toggle1 = Toggle(fig)
    toggle2 = Toggle(fig)
    toggle3 = Toggle(fig)
    toggle4 = Toggle(fig)
    toggle5 = Toggle(fig)

    label1 = Label(fig, "Vorts Unclass")
    label2 = Label(fig, "Vorts Class")
    label3 = Label(fig, "Filter OOB Vorts")
    label4 = Label(fig, "Show Iso")
    label5 = Label(fig, "Vorts Class Lines")

    sl = SliderGrid(fig[2, 1], (label = "Time", range = t[1]:0.001:t[end], format = "{:.1f}s", startvalue = t[1]))

    
    menu_interp = Menu(
        fig, 
        options = collect(1:16),
        default = 1
        )
    menu_label = Label(fig, "Interpolation Depth", width = nothing)

    ϵ = lift(menu_interp.selection) do n
        round(epsilonBall2(sim_vars.X, n, 0.); sigdigits=4)
    end

    
    # epsilon_label = Label(fig, @lift("ε = $(round((ϵ[]); sigdigits=3))"))
    epsilon_label = Label(fig, @lift("ε = $($ϵ) ξ"))
    menu_grid = GridLayout()
    menu_grid[1, 1] = vgrid!(menu_label, menu_interp, epsilon_label, tellheight = false)
    menu_grid[2, 1] = grid!(hcat([toggle1, toggle2, toggle3, toggle4, toggle5], [label1, label2, label3, label4, label5]), tellheight = false)

    # menu_grid[3, 1] = epsilon_label
    # rowgap!(menu_grid, 1, Relative(-0.5))
    fig[1, 2] = menu_grid

    return VortexPlot(fig, [toggle1, toggle2, toggle3, toggle4, toggle5], [label1, label2, label3, label4, label5], [sl.sliders[1]], menu_interp.selection)
end


function density2(psi)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax
    return density
end

mutable struct VortexPlotObservables
    psi::Observable{Array{ComplexF64, 3}}
    density::Observable{Array{Float64, 3}}
    vorts_unclass::Observable
    vorts_class::Observable
    vorts_class_colors::Observable
    vorts_class_lines::Observable
end

function VortexPlotObservables(sim_vars, vortex_plot)
    t = sim_vars.t
    sol = sim_vars.sol
    sl = vortex_plot.sliders[1]
    interp_depth = vortex_plot.interp_depth

    psi = lift(sl.value) do t
        sim_vars.sol(t)
    end

    density = lift(psi) do psi1
        density2(psi1)
    end

    vorts_unclass = VortsUnclassObservable(sim_vars, vortex_plot, psi, interp_depth)

    vorts_class, vorts_class_colors = VortsClassObservable(sim_vars, vortex_plot, psi, interp_depth)

    vorts_class_lines = VortsClassLinesObservable(sim_vars, vortex_plot, psi, interp_depth)
    return VortexPlotObservables(psi, density, vorts_unclass, vorts_class, vorts_class_colors, vorts_class_lines)
end

function PlotVortexPlotObservables(sim_vars, vortex_plot, vortex_plot_observables)
    fig = vortex_plot.fig

    cbarPal= :plasma;
    cmap = get(colorschemes[cbarPal], LinRange(0,1,100));
    cmap2 = [(cmap[i], 0.5) for i in 1:100];

    X = sim_vars.X
    x = X[1]; y = X[2]; z = X[3];

    density = vortex_plot_observables.density
    vorts_unclass = vortex_plot_observables.vorts_unclass
    vorts_class = vortex_plot_observables.vorts_class
    vorts_class_colors = vortex_plot_observables.vorts_class_colors

    show_iso = vortex_plot.toggles[4].active

    volume!(
        fig[1,1], 
        (x[1], x[end]), 
        (y[1], y[end]), 
        (z[1], z[end]), 
        density, 
        algorithm = :iso, 
        isovalue=0.65,
        isorange=0.075, 
        colormap=cmap2, 
        transparency=true,
        visible=show_iso
        )

    meshscatter!(fig[1, 1], vorts_unclass, color=:blue, markersize=0.05)
    meshscatter!(fig[1, 1], vorts_class, color=vorts_class_colors, markersize=0.05)


end

function vortex_points(active, psi, X, N_interp)
    if active
        vorts_3d = find_vortex_points_3d(psi, X, N_interp)
        if length(vorts_3d) == 0
            return Point3f.([], [], [])
        end
        vorts = vorts3DMatrix(vorts_3d)
        return Point3f.(vorts[:, 1], vorts[:, 2], vorts[:, 3])
    else
        return Point3f.([], [], [])
    end
end

function VortsUnclassObservable(sim_vars, vortex_plot, psi, interp_depth)

    vorts_unclass_toggle = vortex_plot.toggles[1]
    vorts_oob_toggle = vortex_plot.toggles[3]
    slider = vortex_plot.sliders[1]
    X = sim_vars.X

    v_points = lift(vorts_unclass_toggle.active, slider.value, interp_depth) do tog, sl, intr
        vortex_points(tog, psi[], X, intr)
    end

    v_points_filt = lift(vorts_oob_toggle.active, v_points) do tog, v
        if tog
            filter(a->vortInBounds(a, X), v)
        else
            v
        end
    end

    return v_points_filt
end

function filter_classified_v(v, X)
    if length(v) == 0
        return Point3f.([], [], [])
    end
    for i in eachindex(v)
        v[i] = filter(a->vortInBounds(a, X), v[i])
    end
    return v
end

function v_color_class(vorts_class, vorts_3d)
    if length(vorts_3d) == 0
        return [RGB(0.8103465454545455,0.2951044545454546,0.4575856363636363)]
    end
    colors = distinguishable_colors(length(vorts_class),[RGB(0.8103465454545455,0.2951044545454546,0.4575856363636363)],dropseed=true);
    v_color = [colors[i] for i in eachindex(vorts_class) for j in eachindex(vorts_class[i])]
    return v_color
end

function arrange_vort_class_for_plotting(vorts_3d, vorts_class)
    if length(vorts_3d) == 0
        return Point3f.([], [], [])
    end
    v_plotting = vcat(vorts_class...)
    v_plotting = reduce(hcat, v_plotting)'
    v_plotting = Point3f.(v_plotting[:, 1], v_plotting[:, 2], v_plotting[:, 3])
    return v_plotting
end

function VortsClassObservable(sim_vars, vortex_plot, psi, interp_depth)
    vorts_class_toggle = vortex_plot.toggles[2]
    vorts_oob_toggle = vortex_plot.toggles[3]
    slider = vortex_plot.sliders[1]
    X = sim_vars.X


    v_points_for_class = lift(vorts_class_toggle.active, slider.value, interp_depth, vorts_oob_toggle.active) do tog, sl, intr, tog2
        if tog
            find_vortex_points_3d(psi[], X, intr)
        else
            []
        end
    end

    v_points_class_2 = lift(v_points_for_class) do v
        connect_vortex_points_3d(v, X, interp_depth[], epsilonBall2(X, interp_depth[], 0.1),  true)
    end

    v_points_class_2_array = lift(v_points_class_2) do v
        [v_points_for_class[][collect(v[i])] for i in eachindex(v)]
    end

    v_points_class_2_array_filt = lift(vorts_oob_toggle.active, v_points_class_2_array) do tog, v
        if tog
            filter_classified_v(v, X)
        else
            v
        end
    end

    v_color = lift(v_points_class_2_array_filt) do v
        v_color_class(v, v_points_for_class[])
    end

    v_plotting = lift(v_points_class_2_array_filt) do v
        arrange_vort_class_for_plotting(v_points_for_class[], v)
    end

    return v_plotting, v_color
end

function VortsClassLinesObservable(sim_vars, vortex_plot, psi, interp_depth)
    vorts_class_lines_toggle = vortex_plot.toggles[5]
    slider = vortex_plot.sliders[1]
    X = sim_vars.X

    v_points_for_class = lift(vorts_class_lines_toggle.active, slider.value, interp_depth) do tog, sl, intr
        if tog
            find_vortex_points_3d(psi[], X, intr)
        else
            []
        end
    end

    v_points_class_2 = lift(v_points_for_class) do v
        connect_vortex_points_3d(v, X, interp_depth[], epsilonBall2(X, interp_depth[], 0.1),  true)
    end

    v_sort_lines = lift(v_points_class_2) do v
        if length(v) == 0
            []
        else
            sort_classified_vorts_3d(v, v_points_for_class[], X)
        end
    end

    v_sort_lines_plotting = lift(v_sort_lines) do v
        if length(v) == 0
            []
        else
            arrange_vort_lines_for_plotting(v, X)
        end
    end

    return v_sort_lines_plotting
end

function arrange_vort_lines_for_plotting(vorts_sorted, X)
    x = X[1]; dx = x[2]-x[1];

    vorts_plotting = []

    for i in 1:length(vorts_sorted)
        vort_current = []
        vi = vorts_sorted[i]
        prev_break = 1
        for j in 1:length(vi)-1
            if euclid(vi[j], vi[j+1]) > dx
                push!(vort_current, vi[prev_break:j])
                prev_break = j+1
            end
        end
        push!(vort_current, vi[prev_break:end])
        if euclid(vort_current[1][1], vort_current[end][end]) < dx
            push!(vort_current[end], vort_current[1][1])
        end

        push!(vorts_plotting, vort_current)
    end
    return vorts_plotting
end

function findVorticesJumps3D(psi, sim)

    X = sim.X
    x = X[1]; y = X[2]; z = X[3];
    x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];

    vorts_xslice = []
    vorts_yslice = []
    vorts_zslice = []

    for i in eachindex(x)
        ψ = Torus(psi[i, :, :], y, z)
        # vorts_x = remove_vortices_edge(findvortices_jumps(ψ), ψ)
        vorts_x = findvortices_jumps(ψ)
        v_slice_temp = []
        for v in vorts_x
            vx = [x[i], v.xv, v.yv]
            push!(v_slice_temp, vx)
        end
        push!(vorts_xslice, v_slice_temp)
    end

    for i in eachindex(y)
        ψ = Torus(psi[:, i, :], x, z)
        # vorts_y = remove_vortices_edge(findvortices_jumps(ψ), ψ)
        vorts_y = findvortices_jumps(ψ)
        v_slice_temp = []
        for v in vorts_y
            vy = [v.xv, y[i], v.yv]
            push!(v_slice_temp, vy)
        end
        push!(vorts_yslice, v_slice_temp)
    end

    for i in eachindex(z)
        ψ = Torus(psi[:, :, i], x, y)
        # vorts_z = remove_vortices_edge(findvortices_jumps(ψ), ψ)
        vorts_z = findvortices_jumps(ψ)
        v_slice_temp = []
        for v in vorts_z
            vz = [v.xv, v.yv, z[i]]
            push!(v_slice_temp, vz)
        end
        push!(vorts_zslice, v_slice_temp)
    end
    return vorts_xslice, vorts_yslice, vorts_zslice
end

function vortex_adjacency_matrix(vorts)
    adj_mat = zeros(Bool, length(vorts), length(vorts))

    for i in eachindex(vorts)
        vi = vorts[i]
        is_x_vort = vi[1] % 0.25 ≈ 0
        is_y_vort = vi[2] % 0.25 ≈ 0
        is_z_vort = vi[3] % 0.25 ≈ 0

        x_lower_bound = is_x_vort ? vi[1] - 0.25 : vi[1] - 0.125
        x_upper_bound = is_x_vort ? vi[1] + 0.25 : vi[1] + 0.125
        y_lower_bound = is_y_vort ? vi[2] - 0.25 : vi[2] - 0.125
        y_upper_bound = is_y_vort ? vi[2] + 0.25 : vi[2] + 0.125
        z_lower_bound = is_z_vort ? vi[3] - 0.25 : vi[3] - 0.125
        z_upper_bound = is_z_vort ? vi[3] + 0.25 : vi[3] + 0.125

        vorts_in_range = findall(x -> 
                        x[1] >= x_lower_bound 
                        && x[1] <= x_upper_bound 
                        && x[2] >= y_lower_bound 
                        && x[2] <= y_upper_bound 
                        && x[3] >= z_lower_bound 
                        && x[3] <= z_upper_bound, 
                        vorts)
                        
        for j in vorts_in_range
            if i != j
                adj_mat[i, j] = true
            end
        end
    end
    return adj_mat
end

function vortex_graph(vorts, adj_mat, graph_plot=true, iso=true)
    g = SimpleGraph(adj_mat)
    if graph_plot
        v_pos = (adj) -> vorts
        if iso
            scene = plot_iso(psi, sim.X)
            graphplot!(g, layout=v_pos, markersize=0.02)
        else
            scene = graphplot(g, layout=v_pos, markersize=0.02)
        end
    else
        scene = Nullscene()
    end
    return g, scene
end

function link_graph_vorts(g; repeat_end_of_ring=false)
    g_temp = deepcopy(g)
    visited = Set()
    deg_1 = Set(findall(x -> x == 1, degree(g_temp)))

    # Link vortex lines by starting at degree 1 vortices (end points) and terminating at degree != 2 vertices
    vort_lines = []
    while length(deg_1) > 0
        current_vort = []
        vc = pop!(deg_1)
        push!(current_vort, vc)
        vc_neighbors = neighbors(g_temp, vc)
        # setdiff!(vc_neighbors, visited)
        while length(vc_neighbors) <= 2 # vc isn't a reconnection
            push!(visited, vc)
            setdiff!(vc_neighbors, visited)
            # print(length(vc_neighbors))
            if length(vc_neighbors) == 0
                break
            end
            vc = pop!(vc_neighbors)
            push!(current_vort, vc)
            vc_neighbors = neighbors(g_temp, vc)  
        end
        push!(vort_lines, current_vort)
        setdiff!(deg_1, visited)
    end

    # Only rings are left, now link rings that have reconnection points
    g_temp = deepcopy(g)
    deg_2g = Set(findall(x -> x > 2, degree(g_temp)))
    vort_loops = []

    while length(deg_2g) > 0
        current_vort = []
        vi = pop!(deg_2g)
        push!(current_vort, vi)
        push!(visited, vi)
        vi_neighbors = neighbors(g_temp, vi)

        # Now traverse a neighbor that is not a reconnection point and not visited
        vi = nothing
        for v in vi_neighbors
            if v ∉ visited
                if length(neighbors(g_temp, v)) <= 2 # This should only be 2 but just in case 
                    vi = v
                    break
                end
            end
        end

        while !(isnothing(vi))
            push!(current_vort, vi)
            push!(visited, vi)
            vi_neighbors = neighbors(g_temp, vi)
            setdiff!(vi_neighbors, visited)
            if length(vi_neighbors) == 1
                vi = pop!(vi_neighbors)
                # Check not reconnection point
                if length(neighbors(g_temp, vi)) > 2
                    push!(current_vort, vi)
                    push!(visited, vi)
                    vi = nothing
                end
            else
                vi = nothing
            end
        end

        # Now we've traversed the ring, remove visited vertices
        push!(vort_loops, current_vort)
        setdiff!(deg_2g, visited)
    end
        
    # Now what's left are rings that don't have reconnection points
    g_temp = deepcopy(g)

    deg_2 = Set(findall(x -> x == 2, degree(g_temp)))
    vort_rings = []
    while length(deg_2) > 0
        current_vort = []
        vi = pop!(deg_2)
        push!(current_vort, vi)
        push!(visited, vi)
        vi_neighbors = neighbors(g_temp, vi)
        setdiff!(vi_neighbors, visited)
        while length(vi_neighbors) != 0
            vc = pop!(vi_neighbors)
            push!(current_vort, vc)
            push!(visited, vc)
            vi_neighbors = neighbors(g_temp, vc)
            setdiff!(vi_neighbors, visited)
        end
        push!(vort_rings, current_vort)
        # put start point at end for ring plotting
        if repeat_end_of_ring
            push!(current_vort, current_vort[1])
        end
        deg_2 = setdiff!(deg_2, visited)
    end

    return vort_lines, vort_loops, vort_rings
    # return vort_lines, vort_loops
end

function smooth_vortex_rings(vort_rings, max_itr=1)
    vort_rings_prev = deepcopy(vort_rings)
    vort_rings_next = deepcopy(vort_rings_prev)

    for itr in 1:max_itr
        for i in eachindex(vort_rings_prev)
            v_ring_length = length(vort_rings_prev[i])
            for j in eachindex(vort_rings_prev[i])
                vort_rings_next[i][j] = ((0.8/5) * sum(vort_rings_prev[i][mod1.(j-2:j+2, v_ring_length)])) + (0.2 * vort_rings_prev[i][j])
                # vort_rings_coords_smoothed_next[i][j] = ((0.8/3) * sum(vort_rings_coords_smoothed_prev[i][mod1.(j-1:j+1, v_ring_length)])) + (0.2 * vort_rings_coords_smoothed_prev[i][j]) # 3 point window
            end
        end
        vort_rings_prev = deepcopy(vort_rings_next)
    end

    return vort_rings_next
end

# Incorporate epsilon tball to then search for grid points

function zoom_vortex_points(vorts_x_slices, vorts_y_slices, vorts_z_slices, psi, sim; periodic=false, win=1)
    X = sim.X
    x = X[1]; y = X[2]; z = X[3];
    x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);

    vorts_x_zoom = []
    for i in eachindex(vorts_x_slices)
        vorts = vorts_x_slices[i]
        ψ = psi[i, :, :]
        
        for (j, vortex) in enumerate(vorts)
            # println(vortex)
            v = 
                try
                psi_int,xint,yint = zoom_interp(ψ, y, z, vortex[2], vortex[3], periodic=periodic, win=win) #periodic: peridic indices here
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))[1]
                psi_int,xint,yint = zoom_interp(psi_int, xint, yint, vint.xv, vint.yv, periodic=periodic, win=win)
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))[1]
                vint = [vortex[1], vint.xv, vint.yv]
                catch
                    nothing
                end
            println("$vortex -> $v")
            isnothing(v) ? push!(vorts_x_zoom, vortex) : push!(vorts_x_zoom, v)
        end
    end


    vorts_y_zoom = []
    for i in eachindex(vorts_y_slices)
        vorts = vorts_y_slices[i]
        ψ = psi[:, i, :]
        
        for (j, vortex) in enumerate(vorts)
            v = try
                psi_int,xint,yint = zoom_interp(ψ, x, z, vortex[1], vortex[3], win=win) #periodic: peridic indices here
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))[1]
                psi_int,xint,yint = zoom_interp(psi_int, xint, yint, vint.xv, vint.yv, win=win)
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))[1]
                vint = [vint.xv, vortex[2], vint.yv]
                catch 
                    nothing
                end
            println("$vortex -> $v")

            isnothing(v) ? push!(vorts_y_zoom, vortex) : push!(vorts_y_zoom, v)
        end
    end

    vorts_z_zoom = []
    for i in eachindex(vorts_z_slices)
        vorts = vorts_z_slices[i]
        ψ = psi[:, :, i]
        
        for (j, vortex) in enumerate(vorts)
            v = try
                psi_int,xint,yint = zoom_interp(ψ, x, y, vortex[1], vortex[2], win=win) #periodic: peridic indices here
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))[1]
                psi_int,xint,yint = zoom_interp(psi_int, xint, yint, vint.xv, vint.yv, win=win)
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))[1]
                vint = [vint.xv, vint.yv, vortex[3]]
                catch 
                    nothing
                end
            println("$vortex -> $v")

            isnothing(v) ? push!(vorts_z_zoom, vortex) : push!(vorts_z_zoom, v)
        end
    end

    return vorts_x_zoom, vorts_y_zoom, vorts_z_zoom
end

function zoom_vortex_points_2(vorts_x_slices, vorts_y_slices, vorts_z_slices, psi, sim; periodic=false, win=0, nz=30)
    X = sim.X
    x = X[1]; y = X[2]; z = X[3];
    x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);

    failed_zoom_count = 0
    vorts_x_zoom = []
    for i in eachindex(vorts_x_slices)
        vorts = vorts_x_slices[i]
        ψ = psi[i, :, :]
        
        for (j, vortex) in enumerate(vorts)
            # println(vortex)
            v = 
                try
                psi_int,xint,yint = zoom_interp(ψ, y, z, vortex[2], vortex[3], periodic=periodic, win=win, nz=nz) #periodic: peridic indices here
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                filter!(a-> a.xv > vortex[2] - 0.125 && a.xv < vortex[2] + 0.125 && a.yv > vortex[3] - 0.125 && a.yv < vortex[3] + 0.125, v1)
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))
                if length(vint) > 0
                    vint = vint[1]
                elseif length(v1) > 0
                    vint = v1[1]
                end
                vint = [vortex[1], vint.xv, vint.yv]
                catch
                    nothing
                end
            # println("$vortex -> $v")
            isnothing(v) ? push!(vorts_x_zoom, vortex) : push!(vorts_x_zoom, v)
            if isnothing(v)
                println("Failed zoom: $vortex")
                failed_zoom_count += 1
            end
        end
    end


    vorts_y_zoom = []
    for i in eachindex(vorts_y_slices)
        vorts = vorts_y_slices[i]
        ψ = psi[:, i, :]
        
        for (j, vortex) in enumerate(vorts)
            v = try
                psi_int,xint,yint = zoom_interp(ψ, x, z, vortex[1], vortex[3], win=win, nz=nz) #periodic: peridic indices here
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                filter!(a-> a.xv > vortex[1] - 0.125 && a.xv < vortex[1] + 0.125 && a.yv > vortex[3] - 0.125 && a.yv < vortex[3] + 0.125, v1)
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))
                if length(vint) > 0
                    vint = vint[1]
                elseif length(v1) > 0
                    vint = v1[1]
                end
                vint = [vint.xv, vortex[2], vint.yv]
                catch 
                    nothing
                end
            # println("$vortex -> $v")

            isnothing(v) ? push!(vorts_y_zoom, vortex) : push!(vorts_y_zoom, v)
            if isnothing(v)
                println("Failed zoom: $vortex")
                failed_zoom_count += 1
            end
        end
    end

    vorts_z_zoom = []
    for i in eachindex(vorts_z_slices)
        vorts = vorts_z_slices[i]
        ψ = psi[:, :, i]
        
        for (j, vortex) in enumerate(vorts)
            v = try
                psi_int,xint,yint = zoom_interp(ψ, x, y, vortex[1], vortex[2], win=win, nz=nz) #periodic: peridic indices here
                v1 = findvortices_grid(Torus(psi_int, xint, yint))
                filter!(a-> a.xv > vortex[1] - 0.125 && a.xv < vortex[1] + 0.125 && a.yv > vortex[2] - 0.125 && a.yv < vortex[2] + 0.125, v1)
                vint = remove_vortices_edge(v1,Torus(psi_int, xint, yint))
                if length(vint) >0
                    vint = vint[1]
                elseif length(v1) > 0
                    vint = v1[1]
                end
                vint = [vint.xv, vint.yv, vortex[3]]
                catch 
                    nothing
                end
            # println("$vortex -> $v")
            isnothing(v) ? push!(vorts_z_zoom, vortex) : push!(vorts_z_zoom, v)
            if isnothing(v)
                println("Failed zoom: $vortex")
                failed_zoom_count += 1
            end
        end
    end
    println("Failed zoom count: $failed_zoom_count")
    return vorts_x_zoom, vorts_y_zoom, vorts_z_zoom
end

