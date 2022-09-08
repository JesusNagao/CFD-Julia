
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

#Método para cambiar el perfil inicial de la ola
function Initwave!(w::Wave,t::String)

    if t == "Gaussian"
        print("Escribe σ y μ en ese orden")
        σ::Float64 = parse(Float64, readline())
        μ::Float64 = parse(Float64, readline())
        w.u[:,1] = Gaussian(w.u,σ,μ)
    end
    
    if t == "Hat"
        w.u = ones(w.nx,w.nt)
    end

end

#Evolución temporal de la ola con la ecuación de difusión
function Difussion(w::Wave; exportcsv = false)
    nx::Int64 = w.nx
    nt::Int64 = w.nt
    ν::Float64 = 0.3
    σ::Float64 = 0.2
    dx::Float64 = 2/(nx-1)
    dt::Float64 = σ*dx^2/ν

    w.u[floor(Int64,0.5/dx):floor(Int64,1/ dx + 1),1] .= 2

    difusse = @animate for i in range(1,stop=nt-1)

        @threads for j in range(2,stop = nx-1)
            w.u[j,i+1] = w.u[j,i] + ν*dt/dx^2*(w.u[j+1,i] - 2*w.u[j,i] + w.u[j-1,i])    
        end
        
        plot(LinRange(0,2,nx),w.u[:,i], xlabel = "u", ylabel = "x", title = "Difusión de una onda", xlims = (0,2), ylims = (0.5,3))

    end
    
    if exportcsv == true
        writedlm("Diff.csv", w.u, ',')
    end

    gif(difusse, "dif.gif", fps = 1000)

end

dif = Wave(100,200)
Initwave!(dif,"Hat")
@time Difussion(dif,exportcsv = true) 