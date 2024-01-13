using GLMakie
using Base.Threads

nx = 41
ny = 41
nt = 500
nit = 50
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = LinRange(0, 2, nx)
y = LinRange(0, 2, ny)

rho = 1.
nu = .1
dt = .001

u = zeros((ny, nx))
v = zeros((ny, nx))
p = zeros((ny, nx)) 
b = zeros((ny, nx))

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

X,Y = meshgrid(x,y)

function build_up_b(b::Array{Float64,2}, rho::Float64, dt::Float64, u::Array{Float64,2}, v::Array{Float64,2}, dx::Float64, dy::Float64)
    
    b[2:end-1, 2:end-1] = (rho * (1 / dt * 
                    ((u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / 
                     (2 * dx) + (v[3:end,2:end-1] - v[1:end-2, 2:end-1]) / (2 * dy)) -
                    ((u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / (2 * dx))^2 -
                      2 * ((u[3:end, 2:end-1] - u[1:end-2, 2:end-1]) / (2 * dy) *
                           (v[2:end-1, 3:end] - v[2:end-1, 1:end-2]) / (2 * dx))-
                          ((v[3:end, 2:end-1] - v[1:end-2, 2:end-1]) / (2 * dy))^2))

    return b
end

function pressure_poisson(p::Array{Float64,2}, dx::Float64, dy::Float64, b::Array{Float64,2})

    pn = copy(p)
    
    @threads for q in range(1, stop = nit)

        pn = copy(p)

        p[2:end-1, 2:end-1] = (((pn[2:end-1, 3:end] + pn[2:end-1, 1:end-2]) * dy^2 + 
                          (pn[3:end, 2:end-1] + pn[1:end-2, 2:end-1]) * dx^2) /
                          (2 * (dx^2 + dy^2)) -
                          dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * 
                          b[2:end-1,2:end-1])

        p[:, end] = p[:, end-1] # dp/dx = 0 at x = 2
        p[1, :] = p[2, :]   # dp/dy = 0 at y = 0
        p[:, 1] = p[:, 2]   # dp/dx = 0 at x = 0
        p[end, :] .= 0        # p = 0 at y = 2
        
    end

    return p

end

function cavity_flow(nt::Int64, u::Array{Float64,2}, v::Array{Float64,2}, dt::Float64, dx::Float64, dy::Float64, p::Array{Float64,2}, rho::Float64, nu::Float64)

    b = zeros((ny, nx))
    
    flow = @animate for n in range(1, stop = nt)
        un = copy(u)
        vn = copy(v)
        
        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b)
        
        u[2:end-1, 2:end-1] = (un[2:end-1, 2:end-1]-
                         un[2:end-1, 2:end-1] * dt / dx *
                        (un[2:end-1, 2:end-1] - un[2:end-1, 1:end-2]) -
                         vn[2:end-1, 2:end-1] * dt / dy *
                        (un[2:end-1, 2:end-1] - un[1:end-2, 2:end-1]) -
                         dt / (2 * rho * dx) * (p[2:end-1, 3:end] - p[2:end-1, 1:end-2]) +
                         nu * (dt / dx^2 *
                        (un[2:end-1, 3:end] - 2 * un[2:end-1, 2:end-1] + un[2:end-1, 1:end-2]) +
                         dt / dy^2 *
                        (un[3:end, 2:end-1] - 2 * un[2:end-1, 2:end-1] + un[1:end-2, 2:end-1])))

        v[2:end-1,2:end-1] = (vn[2:end-1, 2:end-1] -
                        un[2:end-1, 2:end-1] * dt / dx *
                       (vn[2:end-1, 2:end-1] - vn[2:end-1, 1:end-2]) -
                        vn[2:end-1, 2:end-1] * dt / dy *
                       (vn[2:end-1, 2:end-1] - vn[1:end-2, 2:end-1]) -
                        dt / (2 * rho * dy) * (p[3:end, 2:end-1] - p[1:end-2, 2:end-1]) +
                        nu * (dt / dx^2 *
                       (vn[2:end-1, 3:end] - 2 * vn[2:end-1, 2:end-1] + vn[2:end-1, 1:end-2]) +
                        dt / dy^2 *
                       (vn[3:end, 2:end-1] - 2 * vn[2:end-1, 2:end-1] + vn[1:end-2, 2:end-1])))

        u[1, :]  .= 0
        u[:, 1]  .= 0
        u[:, end] .= 0
        u[end, :] .= 1    # set velocity on cavity lid equal to 1
        v[1, :]  .= 0
        v[end, :] .= 0
        v[:,1]  .= 0
        v[:, end] .= 0

        ##How do I plot the resulting X,Y and p fields using GLMakie



        Plots.contourf(x,y,p, xlims = (0,2), ylims = (0,2), c = cgrad(:ice))
        Plots.quiver!(Y,X,quiver=(u[:,:],v[:,:]), xlims = (0,2), ylims = (0,2),c=:viridis)
        
    end

    gif(flow, "2DCavityFlow.gif", fps = 1000)  

    return u, v, p

end

u = zeros((ny, nx))
v = zeros((ny, nx))
p = zeros((ny, nx))
b = zeros((ny, nx))
nt = 700
@time u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu)
