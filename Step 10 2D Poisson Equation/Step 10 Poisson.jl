
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

function Poisson(P::PressureField2D, exportocsv = false)
    nx::Int64 = P.nx
    ny::Int64 = P.ny
    nt::Int64 = 10
    dx::Float64 = 2/(nx-1)
    dy::Float64 = 1/(ny-1) 

    X = LinRange(0,2,nx)
    Y = LinRange(0,1,ny)

    #Condiciones de frontera
    P.p[1,1:end] .= 0
    P.p[end,1:end] .= 0
    P.p[1:end,1] .= 0
    P.p[1:end,end] .= 0    
    
    #Vector con el término no homogeneo
    b = zeros(nx,ny)
    b[ceil(Int64,nx / 4), ceil(Int64,ny / 4)]  = 100
    b[ceil(Int64,3 * nx / 4), ceil(Int64,3 * ny / 4)] = -100
    println(b)
    for i in range(1,stop=nt)
        pcopy = copy(P.p)

        P.p[2:end-1,2:end-1] = ((dy^2 .*(pcopy[2:end-1,3:end] + pcopy[2:end-1,1:end-2]) + dx^2 .* (pcopy[3:end, 2:end-1] + pcopy[1:end-2, 2:end-1]) - (b[2:end-1,2:end-1].*dx^2 .*dy^2))/(2*(dx^2+dy^2)))

        P.p[1,1:end] .= 0
        P.p[end,1:end] .= 0
        P.p[1:end,1] .= 0
        P.p[1:end,end] .= 0    

    end
    
end

function PlotP(P::PressureField2D, exportcsv = false)
    nx::Int64 = P.nx
    ny::Int64 = P.ny

    X = LinRange(0,2,nx)
    Y = LinRange(0,2,ny)

    a = plot(X,Y,P.p[:,:],st = :surface, xlims = (0,2), ylims = (0,2), zlims = (-1,1), camera=(0,20), c = cgrad(:ice))
    savefig(a, "Poiss.png")

    if exportcsv == true
        writedlm("Poiss.csv", P.p, ',')
    end

end

Pos1 = PressureField2D(50,50)
@time Poisson(Pos1)
PlotP(Pos1)