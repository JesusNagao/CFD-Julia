using Plots
using Base.Threads
using DelimitedFiles
using LinearAlgebra

#Paso 12 de Barbara CFD, ecuaciones de momentum de un fluido incompresible

mutable struct NavierStokes2D
    nx::Int64
    ny::Int64
    nt::Int64
    p::Array{Float64,2}   #Presión del fluido
    u::Array{Float64,3}   #Velocidad del elemento de fluido en x
    v::Array{Float64,3}   #Velocidad del elemento de fluido en y
end

NavierStokes2D(nx,ny,nt) = NavierStokes2D(nx,ny,nt,zeros(nx,ny),zeros(nx,ny,nt),zeros(nx,ny,nt)) 

function PoissonP(p::Array{Float64,2},u::Array{Float64,2},v::Array{Float64,2},dx::Float64,dy::Float64,dt::Float64)

    pn::Array{Float64,2} = copy(p)
    nit::Int64 = 50
    ρ::Float64 = 1.

    @threads for i in range(1, stop = nit)
        pn = copy(p)
        
        p[2:end-1, 2:end-1] = ((pn[2:end-1, 3:end] + pn[2:end-1, 1:end-2]) .* dy^2 + 
                              (pn[3:end, 2:end-1] .+ pn[1:end-2, 2:end-1]) .* dx^2) /
                              (2 * (dx^2 + dy^2)) -
                              (ρ * dx^2 * dy^2 / (2 * (dx^2 + dy^2))) .*
                              ((1/ dt) .* ((u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / 
                              (2 * dx) + (v[3:end, 2:end-1] - v[1:end-2, 2:end-1]) / (2 * dy)) -
                              ((u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / (2 * dx)).^2 -
                              2 .* ((u[3:end, 2:end-1] - u[1:end-2, 2:end-1]) / (2 * dy) .*
                              (v[2:end-1, 3:end] - v[2:end-1, 1:end-2]) / (2 * dx)) -
                              ((v[3:end, 2:end-1] - v[1:end-2, 2:end-1]) / (2 * dy)).^2) 
                          
        

        p[1:end,end] = p[1:end, end-2] 
        p[1, 1:end] = p[2, 1:end]   
        p[1:end, 1] = p[1:end, 2]   
        p[end, 1:end] .= 0       

    end

    return p

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

    function meshgrid(x, y)
        X = [i for i in x, j in 1:length(y)]
        Y = [j for i in 1:length(x), j in y]
        return X, Y
    end

    x = LinRange(0,2,nx)
    y = LinRange(0,2,ny)

    X, Y = meshgrid(x,y)

    u::Array{Float64,3} = N.u
    v::Array{Float64,3} = N.v
    p::Array{Float64,2} = N.p

    u[end, 1:end,1] .= 1 

    flow = @animate for i in range(1, stop=nt-1)
        
        p = PoissonP(p,u[:,:,i],v[:,:,i],dx,dy,dt)  #Se asegura la incompresibilidad del campo de densidad cada paso de tiempo
        
        #Resolviendo para cada componente de la ecuación de momentum

        u[2:end-1,2:end-1,i+1] = u[2:end-1, 2:end-1,i] -
                                 u[2:end-1, 2:end-1,i] .* (dt / dx) .*
                                 (u[2:end-1, 2:end-1,i] - u[2:end-1, 1:end-2,i]) -
                                 v[2:end-1, 2:end-1,i] .* (dt / dy) .*
                                 (u[2:end-1, 2:end-1,i] - u[1:end-2, 2:end-1,i]) -
                                 (dt / (2 * ρ * dx)) .* (p[2:end-1, 3:end] - p[2:end-1, 1:end-2]) +
                                 ν .* ((dt / dx^2) .*
                                 (u[2:end-1, 3:end,i] - 2 .* u[2:end-1, 2:end-1,i] + u[2:end-1, 1:end-2,i]) +
                                 (dt / dy^2) .*
                                 (u[3:end, 2:end-1,i] - 2 .* u[2:end-1, 2:end-1,i] + u[1:end-2, 2:end-1,i]))

        v[2:end-1,2:end-1,i+1] = v[2:end-1, 2:end-1,i] -
                                 u[2:end-1, 2:end-1,i] .* (dt / dx) .*
                                 (v[2:end-1, 2:end-1,i] - v[2:end-1, 1:end-2,i]) -
                                 v[2:end-1, 2:end-1,i] .* (dt / dy) .*
                                 (v[2:end-1, 2:end-1,i] - v[1:end-2, 2:end-1,i]) -
                                 (dt / (2 * ρ * dx)) .* (p[2:end-1, 3:end] - p[2:end-1, 1:end-2]) +
                                 ν .* ((dt / dx^2) .*
                                 (v[2:end-1, 3:end,i] - 2 .* v[2:end-1, 2:end-1,i] + v[2:end-1, 1:end-2,i]) +
                                 (dt / dy^2) .*
                                 (v[3:end, 2:end-1,i] - 2 .* v[2:end-1, 2:end-1,i] + v[1:end-2, 2:end-1,i]))
        
        u[1, 1:end,i]  .= 0
        u[1:end, 1,i]  .= 0
        u[1:end, end,i] .= 0
        u[end, 1:end,i] .= 1    
        v[1, 1:end,i]  .= 0
        v[end, 1:end,i] .= 0
        v[1:end, 1,i]  .= 0
        v[1:end, end,i] .= 0

        contourf(x,y,p, xlims = (0,2), ylims = (0,2), c = cgrad(:ice))
        quiver!(X,Y,quiver=(u[:,:,i],v[:,:,i]), xlims = (0,2), ylims = (0,2), aspect_ratio=:equal, arrowscale=0.5)
    end
   
    if exportcsv == true
        writedlm("Navier2Du.csv", u, ',')
        writedlm("Navier2Dv.csv", v, ',')
    end

    gif(flow, "2DCavityFlow.gif", fps = 1000)
    
end

Navier2D = NavierStokes2D(41,41,100)
@time CavityFlow(Navier2D, true)