using CUDA
using GLMakie

CUDA.allowscalar(false)

const nx = 41
const dx = 2 / (nx-1)
const nt = 25
const dt = 0.025
const c = 1

function run()
    u = CuArray{Float64}(undef , nx);
    u = CUDA.ones(nx)

    @btime CUDA.allowscalar() do

        for i in range(Int(round(nx/3)), stop=Int(round(2*nx/3)))

            u[i] = 2;

        end

    end

    up = Array(u) 

    pb = plot(up)

    u_old = CuArray{Float64}(undef, nx)
    u_new = CuArray{Float64}(undef, nx)
    un = CuArray{Float64}(undef, nx-1)

    for i in range(1, stop=nt)
        
        u_old = u[2:nx]
        u_new = u[1:nx-1]

        #println(length(u[2:nx]))

        un = u_old .+ ( - c * dt / dx) * (u_old .- u_new)

        #=
        CUDA.allowscalar() do 

            for i in range(2, stop=nx)
                u[i] = un[i]
            end
        end
        =#
    end

    up = Array(un)

    pa = plot(up)

    plot(pb, pa, layout=(2,1))
end

@btime run()
