function build_up_rho(nx::Int64, ny::Int64)

    for i in range(2, stop=nx-1)
        for j in range(2, stop=ny-1)
            
            b[i+1, j+1] = (rho * (1/dt * ((u[i, j+1]))))

        end
    end



end