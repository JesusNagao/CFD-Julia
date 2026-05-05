
using Plots
using Base.Threads
using DelimitedFiles

#Burger's Equation
#du/dt+u*du/dx=v*d^2u/dx^2

#Create a structure that saves the information of the wave
mutable struct Wave
    nx::Int64
    nt::Int64
    u::Array{Float64,2} 
end

#Instatiate the wave structure
Wave(nx,nt) = Wave(nx,nt,zeros(nx,nt))      

#Define the shape of our wave (Gaussian Shape)
function Gaussian(x::Array{Float64,2},σ::Float64,μ::Float64)
    return [(1/sqrt(2*pi*σ^2)*exp(-((i-μ)^2)/(2*σ^2))) for i in x[:,1]]
end

#Time evolution of the wave using Burger's Equation
function Burgers(w::Wave; exportcsv = false)
    nx::Int64 = w.nx
    nt::Int64 = w.nt
    ν::Float64 = 0.3
    σ::Float64 = 0.2
    dx::Float64 = 2/(nx-1)
    dt::Float64 = σ*dx^2/ν

    #Initialize wave to hat shape
    w.u[floor(Int64,0.5/dx):floor(Int64,1/ dx + 1),1] .= 2

    #Create animation for diffusion for nt steps
    difusse = @animate for i in range(1,stop=nt-1)

        #Run the discretized diffusion equation for the wave
        #@threads is a macro that allows multithreading in julia
        @threads for j in range(2,stop = nx-1)
            w.u[j,i+1] = w.u[j,i] - w.u[j,i]*dt/dx*(w.u[j,i] - w.u[j-1,i]) + ν*dt/dx^2*(w.u[j+1,i] - 2*w.u[j,i] + w.u[j-1,i])    
        end

        #Plot the wave
        plot(LinRange(0,2,nx),w.u[:,i], xlabel = "u", ylabel = "x", title = "Difusión de una onda", xlims = (0,2), ylims = (-3,3))

    end

    #Exporting to CSV
    #if exportcsv == true
    #    writedlm("Burg.csv", w.u, ',')
    #end

    #Create animation
    gif(difusse, "burg.gif", fps = 1000)

end

burg = Wave(100,200)
Initwave!(burg,"Gaussian") #0.2,-1,1
@time Burgers(burg,exportcsv = true) 