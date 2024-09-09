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
    label1 = Label(fig, "Vorts Unclass")
    label2 = Label(fig, "Vorts Class")
    label3 = Label(fig, "Filter OOB Vorts")
    label4 = Label(fig, "Show Iso")

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
    menu_grid[2, 1] = grid!(hcat([toggle1, toggle2, toggle3, toggle4], [label1, label2, label3, label4]), tellheight = false)

    # menu_grid[3, 1] = epsilon_label
    # rowgap!(menu_grid, 1, Relative(-0.5))
    fig[1, 2] = menu_grid

    return VortexPlot(fig, [toggle1, toggle2, toggle3, toggle4], [label1, label2, label3, label4], [sl.sliders[1]], menu_interp.selection)
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

    return VortexPlotObservables(psi, density, vorts_unclass, vorts_class, vorts_class_colors)
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




    