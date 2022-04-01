
using Plots
using Base.Threads
using DelimitedFiles
using LinearAlgebra

mutable struct PressureField2D
    nx::Int64
    ny::Int64
    p::Array{Float64,2} 
end

PressureField2D(nx,ny) = PressureField2D(nx,ny,zeros(nx,ny)) 

function Laplace(P::PressureField2D, exportocsv = false)
    nx::Int64 = P.nx
    ny::Int64 = P.ny
    c::Float64 = 1
    dx::Float64 = 2/(nx-1)
    dy::Float64 = 2/(ny-1) 
    target::Float64 = 0.001
    norml::Float64 = 1.

    #Condiciones de frontera
    P.p[1:end,1] .= 0
    P.p[1:end,end] = Y
    P.p[1,1:end] = P.p[2,1:end]
    P.p[end-1,1:end] .= P.p[end-1,1:end]    
    
    while norml > target
        pcopy = copy(P.p)

        P.p[2:end-1,2:end-1] = ((dy^2 .*(pcopy[2:end-1,3:end] + pcopy[2:end-1,1:end-2]) + dx^2 .* (pcopy[3:end, 2:end-1] + pcopy[1:end-2, 2:end-1]))/(2*(dx^2+dy^2)))

        P.p[1:end,1] .= 0
        P.p[1:end,end] = Y
        P.p[1,1:end] = P.p[2,1:end]
        P.p[end,1:end] .= P.p[end-1,1:end]

        norml = norm(P.p-pcopy)/norm(pcopy)
        println(norml)

    end
end

function PlotP(P::PressureField2D, exportcsv = false)
    nx::Int64 = P.nx
    ny::Int64 = P.ny

    X = LinRange(0,2,nx)
    Y = LinRange(0,2,ny)

    a = plot(X,Y,P.p[:,:],st = :surface, xlims = (0,2), ylims = (0,2), zlims = (0,2), camera=(30,40), c = cgrad(:ice))
    savefig(a, "Lap.png")

    if exportcsv == true
        writedlm("lapl.csv", P.p, ',')
    end

end

Lap1 = PressureField2D(31,31)
@time Laplace(Lap1)
PlotP(Lap1)