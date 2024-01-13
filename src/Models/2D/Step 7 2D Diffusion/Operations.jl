using Plots
#using Threads

function init_speed(u::Matrix{Float64}, x::Array{Float64}, y::Array{Float64})

    for i in range(1, stop=size(u)[1])
        for j in range(1, stop=size(u)[2])
            if (x[i] < 1 && x[i] > 0.5 && y[j] < 1 && y[j] > 0.5 )

                u[i, j] = 2

            end
        end
    end
end

function animate(x::Array{Float64}, y::Array{Float64}, u::Matrix{Float64}, nt::Int64, nx::Int64, ny::Int64, dx::Float64, dy::Float64, dt::Float64, nu::Float64)


    zlims!(0.0, 2.0)

    @gif for k in range(1, stop=nt)
        #s = surface(x, y, u)
        un = u
        diffuse(u, un, nx, ny, dx, dy, dt, nu)
        zlims = (0.0, 2.0)
        surface(x, y, u, xlim=(0,2), ylim=(0,2), zlim=(0,2), legend=false)

    end

end

function diffuse(u::Matrix{Float64}, un::Matrix{Float64}, nx::Int64, ny::Int64, dx::Float64, dy::Float64, dt::Float64, nu::Float64)

    Threads.@threads for i  in range(2, stop=nx-1)
        Threads.@threads for j in range(2, stop=ny-1)
            u[i, j] = un[i,j] + nu*dt/(dx^2)*(un[i+1, j] - 2*un[i,j] + un[i-1, j]) + nu*dt/(dy^2)*(un[i, j+1]-2*un[i,j]+un[i,j-1])
        end
    end

end
