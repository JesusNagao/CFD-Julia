using Plots

function init_speed(u::Matrix{Float64}, x::Array{Float64}, y::Array{Float64})

    for i in range(1, stop=size(u)[1])
        for j in range(1, stop=size(u)[2])
            if (x[i] < 1 && x[i] > 0.5 && y[j] < 1 && y[j] > 0.5 )

                u[i, j] = 2

            end
        end
    end
end

function animate(x::Array{Float64}, y::Array{Float64}, u::Matrix{Float64}, nt::Int64, nx::Int64, ny::Int64)

    s = surface(x, y, u)

    for k in range(1, stop=nt)
        un = u

        for i  in range(1, stop=nx)
            for j in range(1, stop=ny)
                
            end
        end

    end

end
