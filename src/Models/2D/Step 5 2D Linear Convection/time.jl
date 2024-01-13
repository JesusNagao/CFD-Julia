include("space.jl")
using Plots;

function time_iteration(u::Matrix{Float64}, un::Matrix{Float64}, nt::Int64, c::Int64, dt::Float64, dx::Float64, dy::Float64)
    
    @gif for n in range(1,stop=nt+1)
        un = u;
        xlims!(0.0, 2.0)
        ylims!(0.0,2.0)
        zlims!(0.0,2.0)
        surface(x, y, u)
        space_iteration(u, un, c, dt, dx, dy)
    end


end