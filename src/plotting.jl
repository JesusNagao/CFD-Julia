using GLMakie

function plot_3d(x, y, z, title)
    fig = Figure()
    ax = Axis3(fig[1, 1], xlabel = "x", ylabel = "y", zlabel = "z")
    lines!(ax, x, y, z)
    fig[1, 1] = ax
    fig
end

function plot_2d(x, y, title)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y")
    lines!(ax, x, y)
    fig[1, 1] = ax
    fig
end