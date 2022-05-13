using Plots
using Base.Threads
using DelimitedFiles
using LinearAlgebra

#Paso 11 de Barbara CFD, ecuaciones de momentum de un fluido incompresible

mutable struct NavierStokes2D
    nx::Int64
    ny::Int64
    nt::Int64
    p::Array{Float64,3}   #Presión del fluido
    u::Array{Float64,3}   #Velocidad del elemento de fluido en x
    v::Array{Float64,3}   #Velocidad del elemento de fluido en y
end

#Inicializando el objeto
NavierStokes2D(nx,ny,nt) = NavierStokes2D(nx,ny,nt,zeros(nx,ny,nt),zeros(nx,ny,nt),zeros(nx,ny,nt)) 

#Solución de la ecuación de Poisson
function PoissonP(p::Array{Float64,2},u::Array{Float64,2},v::Array{Float64,2},nx::Int64,ny::Int64,dt::Float64)

    pn::Array{Float64,2} = p[:,:]
    nit::Int64 = 100
    dx::Float64 = 2/(nx-1) 
    dy::Float64 = 2/(ny-1) 
    ρ::Float64 = 1.

    @threads for i in range(1, stop = nit)
        pn[:,:] = p[:,:]
        
        p[2:end-1, 2:end-1] = ((pn[2:end-1, 3:end] + pn[2:end-1, 1:end-2]) * dy^2 + 
                              (pn[3:end, 2:end-1] + pn[1:end-2, 2:end-1]) * dx^2) /
                              (2 * (dx^2 + dy^2)) -
                              (ρ * dx^2 * dy^2 / (2 * (dx^2 + dy^2))) *
                              ((1/ dt) * ((u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / 
                              (2 * dx) + (v[3:end, 2:end-1] - v[1:end-2, 2:end-1]) / (2 * dy)) -
                              ((u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / (2 * dx))^2 -
                              2 * ((u[3:end, 2:end-1] - u[1:end-2, 2:end-1]) / (2 * dy) *
                              (v[2:end-1, 3:end] - v[2:end-1, 1:end-2]) / (2 * dx)) -
                              ((v[3:end, 2:end-1] - v[1:end-2, 2:end-1]) / (2 * dy)).^2)
                              
        p[end,:] = p[end-1, :] # dp/dx = 0 at x = 2
        p[:,1] = p[:, 2]       # dp/dy = 0 at y = 0
        p[1, :] = p[2, :]      # dp/dx = 0 at x = 0
        p[:, end] .= 0         # p = 0 at y = 2

        

    end

    return p

end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function CavityFlow(N::NavierStokes2D, exportcsv::Bool)

    nx::Int64 = N.nx
    ny::Int64 = N.ny
    nt::Int64 = N.nt

    c::Float64 = 1.
    ρ::Float64 = 1.
    ν::Float64 = 0.1

    dx::Float64 = 2/(nx-1)
    dy::Float64 = 2/(ny-1) 
    dt::Float64 = 0.0001

    x = LinRange(1,3,nx)
    y = LinRange(1,3,ny)

    X, Y = meshgrid(x,y)

    u::Array{Float64,3} = N.u
    v::Array{Float64,3} = N.v
    p::Array{Float64,3} = N.p

    u[end, :,1] .= 1

    flow = @animate for i in range(1, stop=nt-1)

        
        p[:,:,i] = PoissonP(p[:,:,i],u[:,:,i],v[:,:,i],nx,ny,dt)  #Se asegura la incompresibilidad del campo de densidad cada paso de tiempo  

        #Resolviendo para cada componente de la ecuación de momentum

        u[2:end-1,2:end-1,i+1] = u[2:end-1, 2:end-1,i] -
                                 u[2:end-1, 2:end-1,i] * (dt / dx) *
                                 (u[2:end-1, 2:end-1,i] - u[2:end-1, 1:end-2,i]) -
                                 v[2:end-1, 2:end-1,i] * (dt / dy) *
                                 (u[2:end-1, 2:end-1,i] - u[1:end-2, 2:end-1,i]) -
                                 (dt / (2 * ρ * dx)) * (p[2:end-1, 3:end,i] - p[2:end-1, 1:end-2,i]) +
                                 ν * ((dt / dx^2) *
                                 (u[2:end-1, 3:end,i] - 2 * u[2:end-1, 2:end-1,i] + u[2:end-1, 1:end-2,i]) +
                                 (dt / dy^2) *
                                 (u[3:end, 2:end-1,i] - 2 * u[2:end-1, 2:end-1,i] + u[1:end-2, 2:end-1,i]))

        v[2:end-1,2:end-1,i+1] = v[2:end-1, 2:end-1,i] -
                                 u[2:end-1, 2:end-1,i] * (dt / dx) *
                                 (v[2:end-1, 2:end-1,i] - v[2:end-1, 1:end-2,i]) -
                                 v[2:end-1, 2:end-1,i] * (dt / dy) *
                                 (v[2:end-1, 2:end-1,i] - v[1:end-2, 2:end-1,i]) -
                                 (dt / (2 * ρ * dx)) * (p[2:end-1, 3:end,i] - p[2:end-1, 1:end-2,i]) +
                                 ν * ((dt / dx^2) *
                                 (v[2:end-1, 3:end,i] - 2 * v[2:end-1, 2:end-1,i] + v[2:end-1, 1:end-2,i]) +
                                 (dt / dy^2) *
                                 (v[3:end, 2:end-1,i] - 2 * v[2:end-1, 2:end-1,i] + v[1:end-2, 2:end-1,i]))

        u[1, :,i]  .= 0
        u[:, 1,i]  .= 0
        u[end, :,i] .= 0
        u[:, end,i] .= 1    
        v[1, :,i]  .= 0
        v[:, 1,i]  .= 0
        v[end, :,i] .= 0
        v[:, end,i] .= 0
        


        Plots.contourf(x,y,p[:,:,i], xlims = (0,2), ylims = (0,2), c = cgrad(:turbid))
        Plots.quiver!(X,Y,quiver=(u[:,:,i],v[:,:,i]), xlims = (0,2), ylims = (0,2),c=:viridis)
    end
   
    if exportcsv == true
        writedlm("Navier2Du.csv", u, ',')
        writedlm("Navier2Dv.csv", v, ',')
        writedlm("Navier2Dp.csv", p, ',')
    end

    gif(flow, "2DCavityFlow.gif", fps = 1000)  
    
    return u,v,p
end

Navier2D = NavierStokes2D(41,41,100)
@time u,v,p =  CavityFlow(Navier2D, true)