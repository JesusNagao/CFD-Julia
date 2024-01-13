function init_speed(u::Matrix{Float64}, dx::Float64, dy::Float64)

    for i in range(1, stop=size(u)[1])
        for j in range(1, stop=size(u)[2])

            u[i,j] = 0.0;

            if (j>0.5/dx && j<1/dx+1 && i<1/dy+1 && i>0.5/dy)
                u[i,j] = 2.0;
            end


        end
    end

    

end

function int(x::Float64)

    return floor(Int, x)

end