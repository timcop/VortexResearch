using ColorSchemes, Colors, GLMakie

function plot_iso(psi, X, heal_2=false, visible=true)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax

    cbarPal= :plasma
    cmap = get(colorschemes[cbarPal], LinRange(0,1,100))
    cmap2 = [(cmap[i], 0.5) for i in 1:100]
    xlims = (X[1][1], X[1][end])
    ylims = (X[2][1], X[2][end])
    zlims = (X[3][1], X[3][end])

    if heal_2
        # scene = volume(X[1], X[2], X[3], density, algorithm = :iso, visible=visible, isovalue=0.65,isorange=0.075, colormap=cmap2, transparency=true)
        scene = volume(xlims, ylims, zlims, density, algorithm = :iso, isovalue=0.25,isorange=0.075*2, colormap=cmap2, transparency=true, visible=visible)
    else
        # scene = volume(X[1], X[2], X[3], density, algorithm = :iso, visible=visible, colormap=cmap2, transparency=true, isovalue=0.1,isorange=0.075)
        # scene = volume(density, algorithm = :iso, colormap=cmap2, transparency=true)
        scene = volume(xlims, ylims, zlims, density, algorithm = :iso, colormap=cmap2, transparency=true, visible=visible)
    end

    # scene = volume(X[1], X[2], X[3], density, algorithm = :iso, isovalue=0.9,isorange=0.075, colormap=cmap2, transparency=true)

    # screen = display(scene)


    # resize!(screen, 2998, 1920)
    return scene
end

function vorts3DMatrix(vorts)
    vorts = vcat(vorts'...)
end

function scatterVortsOnIso(ax, vorts, markersize=0.1, transparency=false, alpha=1)
    vorts = vorts3DMatrix(vorts);
    
    meshscatter!(ax, vorts[:, 1], vorts[:, 2], vorts[:, 3], color="blue", markersize=markersize, transparency=transparency, alpha=alpha)
end

function scatterVortsOnIso!(vorts; markersize=0.1, transparency=false, alpha=1, shininess=32.0, color=:blue)
    vorts = vorts3DMatrix(vorts);
    
    meshscatter!(vorts[:, 1], vorts[:, 2], vorts[:, 3], color=color, markersize=markersize, transparency=transparency, alpha=alpha, shininess=shininess, specular=4)
end

function scatterClassifiedVortices(vortSets, vorts_3d, X, edges=false, markersize=0.1)
    colors = distinguishable_colors(length(vortSets),[RGB(0.8103465454545455,0.2951044545454546,0.4575856363636363)],dropseed=true)
    v_matrix = vcat(vorts_3d'...)[:,1:3]'

    for i in 1:length(vortSets)
        print("i")
        vi = v_matrix[:, collect(vortSets[i])]
        if !edges
            vi = vi[:, [vortInBounds(vi[:, i], X) for i = 1:length(vi[1, :])]] # Filters vortices that aren't on the grid
        end
        meshscatter!(vi[1,:],vi[2,:],vi[3,:],color=colors[i], markersize=markersize)
    end
end

# function scatterClassifiedVortices(vortSets, vorts_3d, X, color, edges=false, markersize=0.1)
#     # colors = distinguishable_colors(length(vortSets),[RGB(0.8103465454545455,0.2951044545454546,0.4575856363636363)],dropseed=true)
#     v_matrix = vcat(vorts_3d'...)[:,1:3]'

#     for i in 1:length(vortSets)
#         print("i")
#         vi = v_matrix[:, collect(vortSets[i])]
#         if !edges
#             vi = vi[:, [vortInBounds(vi[:, i], X) for i = 1:length(vi[1, :])]] # Filters vortices that aren't on the grid
#         end
#         meshscatter!(vi[1,:],vi[2,:],vi[3,:],color=color, markersize=markersize)
#     end
# end

function scatterClassifiedVorticesObservable!(vortSets_o, vorts_3d_o, X, edges=false)
    vortSets = vortSets_o[]
    vorts_3d = vorts_3d_o[]
    
    colors = distinguishable_colors(length(vortSets),[RGB(0.8103465454545455,0.2951044545454546,0.4575856363636363)],dropseed=true)
    v_matrix = vcat(vorts_3d'...)[:,1:3]'

    for i in 1:length(vortSets)
        vi = v_matrix[:, collect(vortSets[i])]
        if !edges
            vi = vi[:, [vortInBounds(vi[:, i], X) for i = 1:length(vi[1, :])]] # Filters vortices that aren't on the grid
        end
        meshscatter!(vi[1,:],vi[2,:],vi[3,:],color=colors[i])
    end
end

function vortInBounds(v, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    if ((v[1] >= x[1]) && (v[1] <= x[end]) && 
        (v[2] >= y[1]) && (v[2] <= y[end]) && 
        (v[3] >= z[1]) && (v[3] <= z[end]))
        return true
    else 
        return false
    end
end

function plot_line(vort_sort, color, linewidth)
    vx = [vort_sort[i][1] for i in 1:length(vort_sort)]
    vy = [vort_sort[i][2] for i in 1:length(vort_sort)]
    vz = [vort_sort[i][3] for i in 1:length(vort_sort)]

    lines!(vx, vy, vz, linewidth = linewidth, color = color)
end

function periodicPlotting(vorts_sorted, X, linewidth=5)
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
    
    colors = distinguishable_colors(length(vorts_plotting),[RGB(0.8103465454545455,0.2951044545454546,0.4575856363636363)],dropseed=true);
    for i in 1:length(vorts_plotting)
        vi = vorts_plotting[i]
        try
            for j in 1:length(vi)
                plot_line(vi[j], colors[i], linewidth)
            end
        catch
        end
    end
end

function euclid(v1, v2)
    @assert length(v1) == length(v2)
    sum = 0
    for i in 1:length(v1)
        sum += (v1[i]-v2[i])^2
    end
    sum = sqrt(sum)
    return sum
end