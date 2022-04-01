

function init(u::Matrix{Float64}, v::Matrix{Float64}, nx::Int64, ny::Int64)

    for i in range(1, stop=nx)
        for j in range(1, stop=ny)

            for i in range(1, stop=size(u)[1])
                for j in range(1, stop=size(u)[2])
                    if (x[i] < 1 && x[i] > 0.5 && y[j] < 1 && y[j] > 0.5 )
        
                        u[i, j] = 2
                        v[i, j] = 2
        
                    end
                end
            end

        end
    end

end