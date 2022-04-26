using GLMakie
include("Operations.jl")

const nx = 41
const ny = 41
const nt = 10
const nit = 51
const c = 1
const dx = 2 / (nx - 1)
const dy = 2 / (ny - 1)

x = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = Array{Float64}([i for i in range(0.0, stop=2.0, step=dy)])

const rho = 1
const nu = 0.1
const F = 1
const dt = 0.01

u = Matrix{Float64}(undef, nx, ny)
u = zeros(nx, ny)
v = Matrix{Float64}(undef, nx, ny)
v = zeros(nx, ny)
p = Matrix{Float64}(undef, nx, ny)
p = zeros(nx, ny)
b = Matrix{Float64}(undef, nx, ny)
b = zeros(nx, ny)
un = Matrix{Float64}(undef, nx, ny)
un = zeros(nx, ny)
vn = Matrix{Float64}(undef, nx, ny)
vn = zeros(nx, ny)

X, Y = meshgrid(x, y)

udiff = 1
stepcount = 0

f = Figure(resolution = (800, 800))
Axis(f[1, 1], backgroundcolor = "black")
strength = vec(sqrt.(u .^ 2 .+ v .^ 2))

frames = 1:0.01:10
record(f, "Channel Flow.mp4", frames) do frame
#while (udiff>0.001)
    #if (udiff>0.001)
    un = u
    vn = v

    build_up_b(rho, dt, dx, dy, u, v, b, nx, ny)
    pressure_poisson_periodic(p, dx, dy, nit, nx, ny)

    u[2:nx-1, 2:ny-1] = (un[2:nx-1, 2:ny-1] .- un[2:nx-1, 2:ny-1] * dt / dx * (un[2:nx-1, 2:ny-1] .- un[2:nx-1, 1:ny-2]) .- vn[2:nx-1, 2:ny-1] * dt / dy * (un[2:nx-1, 2:nx-1] .- un[1:nx-2, 2:ny-1]) .- dt / (2 * rho * dx) * (p[2:nx-1, 3:ny] .- p[2:nx-1, 1:ny-2]) .+ nu * (dt / dx^2 * (un[2:nx-1, 3:ny] .- 2 * un[2:nx-1, 2:ny-1] .+ un[2:nx-1, 1:ny-2]) .+ dt / dy^2 * (un[3:nx, 2:ny-1] .- 2 * un[2:nx-1, 2:ny-1] .+ un[1:nx-2, 2:ny-1])) .+ F * dt)
    v[2:nx-1, 2:ny-1] = (vn[2:nx-1, 2:ny-1] .- un[2:nx-1, 2:ny-1] .* dt / dx .* (vn[2:nx-1, 2:ny-1] .- vn[2:nx-1, 1:ny-2]) .- vn[2:nx-1, 2:ny-1] .* dt / dy .* (vn[2:nx-1, 2:ny-1] .- vn[1:nx-2, 2:ny-1]) .- dt / (2 * rho * dy) .* (p[3:nx, 2:ny-1] - p[1:nx-2, 2:ny-1]) .+ nu .* (dt / dx^2 * (vn[2:nx-1, 3:ny] .- 2 .* vn[2:nx-1, 2:ny-1] .+ vn[2:nx-1, 1:ny-2]) .+ dt / dy^2 .* (vn[3:nx, 2:ny-1] .- 2 .* vn[2:nx-1, 2:ny-1] .+ vn[1:nx-2, 2:ny-1])))
    u[2:nx-1, ny-1] = (un[2:nx-1, ny-1] - un[2:nx-1, ny-1] .* dt / dx .* (un[2:nx-1, ny-1] - un[2:nx-1, ny-2]) - vn[2:nx-1, ny-1] .* dt / dy .* (un[2:nx-1, ny-1] - un[1:nx-2, ny-1]) - dt / (2 * rho * dx) .* (p[2:nx-1, 1] - p[2:nx-1, ny-2]) .+ nu .* (dt / dx^2 * (un[2:nx-1, 1] - 2 * un[2:nx-1,ny-1] + un[2:nx-1, ny-2]) + dt / dy^2 .* (un[3:nx, ny-1] - 2 * un[2:nx-1, ny-1] + un[1:nx-2, ny-1])) .+ F .* dt)
    u[2:nx-1, 1] = (un[2:nx-1, 1] - un[2:nx-1, 1] .* dt / dx .* (un[2:nx-1, 1] - un[2:nx-1, ny-1]) - vn[2:nx-1, 1] .* dt / dy .* (un[2:nx-1, 1] - un[1:nx-2, 1]) - dt / (2 * rho * dx) .* (p[2:nx-1, 2] - p[2:nx-1, ny-1]) .+ nu .* (dt / dx^2 * (un[2:nx-1, 2] - 2 * un[2:nx-1, 1] + un[2:nx-1, ny]) + dt / dy^2 .* (un[3:nx, 1] - 2 .* un[2:nx-1, 1] + un[1:nx-2, ny])) .+ F .* dt)
    v[2:nx-1, ny-1] = (vn[2:nx-1, ny-1] - un[2:nx-1, ny-1] .* dt / dx .* (vn[2:nx-1, ny-1] - vn[2:nx-1, ny-2]) - vn[2:nx-1, ny-1] .* dt / dy .* (vn[2:nx-1, ny-1] - vn[1:nx-2, ny-1]) - dt / (2 * rho * dy) .* (p[3:nx, ny-1] - p[1:nx-2, ny-1]) .+ nu .* (dt / dx^2 .* (vn[2:nx-1, 1] - 2 .* vn[2:nx-1, ny-1] + vn[2:nx-1, ny-2]) + dt / dy^2 .* (vn[3:nx, ny-1] - 2 .* vn[2:nx-1, ny-1] + vn[1:nx-2, ny-1])))
    v[2:nx-1, 1] = (vn[2:nx-1, 1] - un[2:nx-1, 1] .* dt / dx .* (vn[2:nx-1, 1] - vn[2:nx-1, ny-1]) - vn[2:nx-1, 1] .* dt / dy .* (vn[2:nx-1, 1] - vn[1:nx-2, 1]) - dt / (2 * rho * dy) .* (p[3:nx, 1] - p[1:nx-2, 1]) .+ nu .* (dt / dx^2 * (vn[2:nx-1, 2] - 2 .* vn[2:nx-1, 1] + vn[2:nx-1, ny-1]) + dt / dy^2 .* (vn[3:nx, 1] - 2 .* vn[2:nx-1, 1] + vn[1:nx-2, 1])))

    u[1, :] .= 0
    u[nx, :] .= 0
    v[1, :] .= 0
    v[nx, :] .= 0
     
    
    #println(stepcount)
    #udiff = (sum(u) - sum(un)) / sum(u)
    #global stepcount += 1
    #end
    #contourf(x,y,p, xlims = (0,2), ylims = (0,2), c = cgrad(:ice))
    #quiver!(X,Y, quiver= u,v , xlims = (0,2), ylims = (0,2), aspect_ratio=:equal, arrowscale=0.5)
    arrows!(x, y, u, v, arrowsize = 10, lengthscale = 0.3,
    arrowcolor = strength, linecolor = strength)
    
end


