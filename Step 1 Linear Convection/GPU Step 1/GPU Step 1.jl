using CUDA
using Plots
include("Operations.jl")
using BenchmarkTools

CUDA.allowscalar(false)

const nx = 41
const dx = 2 / (nx-1)
const nt = 25
const dt = 0.025
const c = 1

function run()

    u = CuArray{Float64}(undef , nx);
    u = CUDA.ones(nx)

    @cuda threads=41 initialize!(u, nx)
    #synchronize()
    
    
    up = Array(u) 
    pb = plot(up)


    u_old = CuArray{Float64}(undef, nx)

    for j in range(1, stop=nt)

        u_old = copy(u)

        @cuda threads=41 kernel_step_1!(u, u_old, c, dt, dx, nx)
        #synchronize()

    end

    up = Array(u)
    pa = plot(up)

    plot(pb, pa, layout=(2,1))

end

@btime run()
#run()