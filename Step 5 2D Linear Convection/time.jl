include("space.jl")
using Plots; pyplot()

function time_iteration(u::Matrix{Float64}, un::Matrix{Float64}, nt::Int64, c::Int64, dt::Float64, dx::Float64, dy::Float64)

    @gif for n in range(1,stop=nt+1)
        space_iteration(u, un, c, dt, dx, dy)
        surface(x, y, u)
    end


end