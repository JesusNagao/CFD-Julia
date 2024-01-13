function space_iteration(u::Matrix{Float64}, un::Matrix{Float64}, c::Int64, dt::Float64, dx::Float64, dy::Float64)
    nx = size(u)[1]
    ny = size(u)[2]

    @threads for i in range(2, stop=nx)
        @threads for j in range(2, stop=ny)
            u[i, j] = (un[i, j] - (c * dt / dx * (un[i, j] - un[i, j-1]))-(c * dt / dy * (un[i, j] - un[i-1, j])))
        end
    end
end