
using Plots
using Base.Threads
using DelimitedFiles
using BenchmarkTools

#Diffusion equation
#du/dt=vd^2u/dt^2

#Create a structure that saves the information of the wave
mutable struct Wave
    nx::Int64
    nt::Int64
    u::Array{Float64,2} 
end

#Instatiate the wave structure
Wave(nx,nt) = Wave(nx,nt,zeros(nx,nt))      

#Define the shape of our wave
function Gaussian(x::Array{Float64,2},σ::Float64,μ::Float64)
    return [(1/sqrt(2*pi*σ^2)*exp(-((i-μ)^2)/(2*σ^2))) for i in x[:,1]]
end

#Initialize wave
function Initwave!(w::Wave,t::String)

    #Use Gaussian distribution as the wave shape
    if t == "Gaussian"
        print("Escribe σ y μ en ese orden")
        σ::Float64 = parse(Float64, readline())
        μ::Float64 = parse(Float64, readline())
        w.u[:,1] = Gaussian(w.u,σ,μ)
    end
    
    #Use Hat distribution as the wave shape
    if t == "Hat"
        w.u = ones(w.nx,w.nt)
    end

end

#Apply diffusion to the wave
function Difussion(w::Wave; exportcsv = false)
    
    #Initialize our variables
    nx::Int64 = w.nx #number of steps in the spatial domain
    nt::Int64 = w.nt #number of time steps
    ν::Float64 = 0.3 #Parameter for our gaussian function
    σ::Float64 = 0.2 #Parameter for our gaussian function
    dx::Float64 = 2/(nx-1) #distance between grid points
    dt::Float64 = σ*dx^2/ν #time step


    w.u[floor(Int64,0.5/dx):floor(Int64,1/ dx + 1),1] .= 2

    #Create animation for diffusion for nt steps
    difusse = @animate for i in range(1,stop=nt-1)

        #Run the discretized diffusion equation for the wave
        #@threads is a macro that allows multithreading in julia
        @threads for j in range(2,stop = nx-1)
            w.u[j,i+1] = w.u[j,i] + ν*dt/dx^2*(w.u[j+1,i] - 2*w.u[j,i] + w.u[j-1,i])    
        end
        
        #Plot the wave
        plot(LinRange(0,2,nx),w.u[:,i], xlabel = "u", ylabel = "x", title = "Difusión de una onda", xlims = (0,2), ylims = (0.5,3))

    end

    #Exporting to CSV
    #if exportcsv == true
    #    writedlm("Diff.csv", w.u, ',')
    #end
    
    
    gif(difusse, "dif.gif", fps = 1000)

end

dif = Wave(100,200)
Initwave!(dif,"Hat")
@btime Difussion(dif,exportcsv = true) 