function init_speed(u::Matrix{Float64}, dx::Float64, dy::Float64)

    for i in range(int(0.5/dy), stop=int(1/dy+1))
        for j in range(int(0.5/dx), stop=int(1/dx+1))
            
            u[i,j] = 2

        end
    end


end

function int(x::Float64)

    return floor(Int, x)


end