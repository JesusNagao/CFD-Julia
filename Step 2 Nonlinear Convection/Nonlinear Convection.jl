using Plots
using BenchmarkTools

#Non-linear Convection Equation
#du/dt+u*du/dx=0

#Initialize variables
const nx = 61 #number of steps in the spatial domain
const dx = 2 / (nx-1) #distance between grid points
const nt = 25 #number of time steps
const dt = 0.025 #time step
const c = 1 #wave propagation speed

function run()

    #Initialize our wave with values of 1
    u = Array{Float64}(undef , nx);
    u = ones(nx)

    #Create a hat pattern
    for i in range(Int(round(nx/3)), stop=Int(round(2*nx/3)))

        u[i] = 2

    end

    #Plot the initial state of the wave
    pb = plot(u, 0, 60)

    #Apply convection for nt steps
    for i in range(1, stop=nt)

        #Copy the previous state of the wave
        u_old = Array{Float64}(undef, nx)
        u_old = copy(u)

        for j in range(2, stop=nx)
            #Calculate the next state of the wave
            u[j] = u_old[j] - u_old[j] * dt / dx * (u_old[j] - u_old[j-1])

        end
        
    end

    #Plot the wave after nt time has passed
    pa = plot(u, 0, 60)
    plot(pb, pa, layout=(2,1))

end

@btime run()