using Plots
using Base.Threads
include("GPUKernel.jl")

function meshgrid(x,y,z)
    X = [i for i in x, j in 1:length(y), k in 1:length(z)]
    Y = [j for i in 1:length(x), j in y, k in 1:length(z)]
    Z = [k for i in 1:length(x), j in 1:length(y), k in z]
    return X, Y, Z
end


function build_up_b(b::Array{Float64,3}, rho::Float64, dt::Float64, u::Array{Float64,3}, v::Array{Float64,3},w::Array{Float64,3}, dx::Float64, dy::Float64, dz::Float64)
    
    b[2:end-1, 2:end-1,2:end-1] = (rho * (1 / dt * 
                    ((u[3:end, 2:end-1,2:end-1] - u[1:end-2, 2:end-1,2:end-1]) / (2 * dx) + 
                     (v[2:end-1,3:end,2:end-1] - v[2:end-1, 1:end-2,2:end-1]) / (2 * dy) +
                     (w[2:end-1,2:end-1,3:end] - w[2:end-1, 2:end-1,1:end-2]) / (2 * dz)) -

                    ((u[3:end, 2:end-1,2:end-1] - u[1:end-2, 2:end-1,2:end-1]) / (2 * dx)).^2 -
                    ((v[2:end-1, 3:end,2:end-1] - v[2:end-1, 1:end-2,2:end-1]) / (2 * dy)).^2 -
                    ((w[2:end-1, 2:end-1,3:end] - w[2:end-1, 2:end-1,1:end-2]) / (2 * dz)).^2 -

                      2 * ((u[2:end-1, 3:end,2:end-1] - u[2:end-1, 1:end-2,2:end-1]) / (2 * dx) .*
                           (v[3:end, 2:end-1,2:end-1] - v[1:end-2, 2:end-1,2:end-1]) / (2 * dy))-
                    
                      2 * ((w[3:end, 2:end-1,2:end-1] - w[1:end-2, 2:end-1,2:end-1]) / (2 * dx) .*
                           (u[2:end-1, 2:end-1, 3:end] - u[2:end-1, 2:end-1,1:end-2]) / (2 * dz))-

                      2 * ((v[2:end-1, 2:end-1,3:end] - v[2:end-1, 2:end-1,1:end-2]) / (2 * dz) .*
                           (w[2:end-1, 3:end,2:end-1] - w[2:end-1, 1:end-2,2:end-1]) / (2 * dy))))

    return b
end

function pressure_poisson(p::Array{Float64,3}, dx::Float64, dy::Float64, dz::Float64, b::Array{Float64,3})

    pn = copy(p)
    
    @threads for q in range(1, stop = nit)

        pn = copy(p)

        p[2:end-1, 2:end-1,2:end-1] = (((pn[3:end, 2:end-1,2:end-1] + pn[1:end-2, 2:end-1,2:end-1]) * dx^2 + 
                                        (pn[2:end-1, 3:end,2:end-1] + pn[2:end-1, 1:end-2,2:end-1]) * dy^2 +
                                        (pn[2:end-1, 2:end-1,3:end] + pn[2:end-1, 2:end-1,1:end-2]) * dz^2) /
                                        (2 * (dx^2 + dy^2 + dz^2)) -
                                        (dx^2 * dy^2 * dz^2) / (2 * (dx^2 + dy^2+ dz^2)) * 
                                        b[2:end-1,2:end-1,2:end-1])

        p[end, :,:] = p[end-1, :,:] # dp/dx = 0 at x = 2
        p[:, 1,:] = p[:, 2,:]   # dp/dy = 0 at y = 0
        p[1, :,:] = p[2, :,:]   # dp/dx = 0 at x = 0
        p[:, end,:] .= 0        # p = 0 at y = 2     
        
    end

    return p

end

function cavity_flow_threads(nt::Int64, u::Array{Float64,3}, v::Array{Float64,3},w::Array{Float64,3}, dt::Float64, dx::Float64, dy::Float64,dz::Float64, p::Array{Float64,3}, rho::Float64, nu::Float64)

    b = CuArray(zeros((nx, ny, nz)))
    
    @threads for n in range(1, stop = nt)

        un = copy(u)
        vn = copy(v)
        wn = copy(w)
        
        b = build_up_b(b, rho, dt, u, v, w, dx, dy, dz)
        p = pressure_poisson(p, dx, dy, dz, b)
        
        u[2:end-1, 2:end-1,2:end-1] = (un[2:end-1, 2:end-1,2:end-1]-
                         un[2:end-1, 2:end-1,2:end-1] * dt / dx .*
                        (un[2:end-1, 2:end-1,2:end-1] - un[1:end-2, 2:end-1,2:end-1]) -
                         vn[2:end-1, 2:end-1,2:end-1] * dt / dy .*
                        (un[2:end-1, 2:end-1,2:end-1] - un[2:end-1, 1:end-2,2:end-1]) -
                         wn[2:end-1, 2:end-1,2:end-1] * dt / dz .*
                        (un[2:end-1, 2:end-1,2:end-1] - un[2:end-1, 2:end-1,1:end-2]) -
                         dt / (2 * rho * dx) * (p[3:end, 2:end-1,2:end-1] - p[1:end-2, 2:end-1,2:end-1]) +
                         nu * 
                         (dt / dx^2 *
                        (un[3:end, 2:end-1,2:end-1] - 2 * un[2:end-1, 2:end-1,2:end-1] + un[1:end-2, 2:end-1,2:end-1]) +
                         dt / dy^2 *
                        (un[2:end-1, 3:end,2:end-1] - 2 * un[2:end-1, 2:end-1,2:end-1] + un[2:end-1, 1:end-2,2:end-1]) +
                        dt / dz^2 *
                        (un[2:end-1, 2:end-1,3:end] - 2 * un[2:end-1, 2:end-1,2:end-1] + un[2:end-1, 2:end-1,1:end-2])))

        v[2:end-1, 2:end-1,2:end-1] = (vn[2:end-1, 2:end-1,2:end-1]-
                         un[2:end-1, 2:end-1,2:end-1] * dt / dx .*
                        (vn[2:end-1, 2:end-1,2:end-1] - vn[1:end-2, 2:end-1,2:end-1]) -
                         vn[2:end-1, 2:end-1,2:end-1] * dt / dy .*
                        (vn[2:end-1, 2:end-1,2:end-1] - vn[2:end-1, 1:end-2,2:end-1]) -
                         wn[2:end-1, 2:end-1,2:end-1] * dt / dz .*
                        (vn[2:end-1, 2:end-1,2:end-1] - vn[2:end-1, 2:end-1,1:end-2]) -
                         dt / (2 * rho * dy) * (p[2:end-1, 3:end,2:end-1] - p[2:end-1, 1:end-2,2:end-1]) +
                         nu * 
                         (dt / dx^2 *
                        (vn[3:end, 2:end-1,2:end-1] - 2 * vn[2:end-1, 2:end-1,2:end-1] + vn[1:end-2, 2:end-1,2:end-1]) +
                         dt / dy^2 *
                        (vn[2:end-1, 3:end,2:end-1] - 2 * vn[2:end-1, 2:end-1,2:end-1] + vn[2:end-1, 1:end-2,2:end-1]) +
                        dt / dz^2 *
                        (vn[2:end-1, 2:end-1,3:end] - 2 * vn[2:end-1, 2:end-1,2:end-1] + vn[2:end-1, 2:end-1,1:end-2])))

        
        w[2:end-1, 2:end-1,2:end-1] = (wn[2:end-1, 2:end-1,2:end-1]-
                         un[2:end-1, 2:end-1,2:end-1] * dt / dx .*
                        (wn[2:end-1, 2:end-1,2:end-1] - wn[1:end-2, 2:end-1,2:end-1]) -
                         vn[2:end-1, 2:end-1,2:end-1] * dt / dy .*
                        (wn[2:end-1, 2:end-1,2:end-1] - wn[2:end-1, 1:end-2,2:end-1]) -
                         wn[2:end-1, 2:end-1,2:end-1] * dt / dz .*
                        (wn[2:end-1, 2:end-1,2:end-1] - wn[2:end-1, 2:end-1,1:end-2]) -
                         dt / (2 * rho * dz) * (p[2:end-1, 2:end-1,3:end] - p[2:end-1, 2:end-1,1:end-2]) +
                         nu * 
                         (dt / dx^2 *
                        (vn[3:end, 2:end-1,2:end-1] - 2 * vn[2:end-1, 2:end-1,2:end-1] + vn[1:end-2, 2:end-1,2:end-1]) +
                         dt / dy^2 *
                        (vn[2:end-1, 3:end,2:end-1] - 2 * vn[2:end-1, 2:end-1,2:end-1] + vn[2:end-1, 1:end-2,2:end-1]) +
                        dt / dz^2 *
                        (vn[2:end-1, 2:end-1,3:end] - 2 * vn[2:end-1, 2:end-1,2:end-1] + vn[2:end-1, 2:end-1,1:end-2])))

        u[1, :,:]  .= 0
        u[:, 1,:]  .= 0
        u[:, end,:] .= 1
        u[end, :,:] .= 0   
        u[:,:,1] .= 0
        u[:,:,end] .= 0

        v[1, :,:]  .= 0
        v[end, :,:] .= 0
        v[:,1,:]  .= 0
        v[:, end,:] .= 0
        v[:,:,1] .= 0
        v[:,:,end] .= 0

        w[1, :,:]  .= 0
        w[end, :,:] .= 0
        w[:,1,:]  .= 0
        w[:, end,:] .= 0
        w[:,:,1] .= 0
        w[:,:,end] .= 0
        
    end

    return u, v, w, p

end
