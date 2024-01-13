
using Plots
using Base.Threads
using DelimitedFiles

mutable struct Wave
    nx::Int64
    nt::Int64
    u::Array{Float64,2} 
end

#Método para inicializar el tipo solo con el número de puntos en x y en t
Wave(nx,nt) = Wave(nx,nt,zeros(nx,nt))      

function Gaussian(x::Array{Float64,2},σ::Float64,μ::Float64)
    return [(1/sqrt(2*pi*σ^2)*exp(-((i-μ)^2)/(2*σ^2))) for i in x[:,1]]
end

#Evolución temporal de la ola con la ecuación de Burgers
function Burgers(w::Wave; exportcsv = false)
    nx::Int64 = w.nx
    nt::Int64 = w.nt
    ν::Float64 = 0.3
    σ::Float64 = 0.2
    dx::Float64 = 2/(nx-1)
    dt::Float64 = σ*dx^2/ν

    w.u[floor(Int64,0.5/dx):floor(Int64,1/ dx + 1),1] .= 2

    difusse = @animate for i in range(1,stop=nt-1)

        @threads for j in range(2,stop = nx-1)
            w.u[j,i+1] = w.u[j,i] - w.u[j,i]*dt/dx*(w.u[j,i] - w.u[j-1,i]) + ν*dt/dx^2*(w.u[j+1,i] - 2*w.u[j,i] + w.u[j-1,i])    
        end
        
        plot(LinRange(0,2,nx),w.u[:,i], xlabel = "u", ylabel = "x", title = "Difusión de una onda", xlims = (0,2), ylims = (-3,3))

    end

    if exportcsv == true
        writedlm("Burg.csv", w.u, ',')
    end

    gif(difusse, "burg.gif", fps = 1000)

end

burg = Wave(100,200)
Initwave!(burg,"Gaussian") #0.2,-1,1
@time Burgers(burg,exportcsv = true) 