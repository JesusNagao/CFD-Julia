function pressure_equation()

    for i in range(1, stop=nx)
        for j in range(1, stop=ny)

            p[i,j] = ((p[i+1, j] + p[i-1, j])*(dy^2) + (p[i, j+1] + p[i, j-1])*(dx^2))/(2*((dx^2)+(dy^2))) - (rho*(dx^2)*(dy^2)/(2*((dx^2) + (dy^2))))*(((u[i+1,j] - u[i-1, j])/2*dx + (v[i,j+1] - v[i, j-1])/2*dy)/dt-((u[i+1, j] - u[i-1, j])*(u[i+1, j]-u[i-1,j])/(4*dx^2))-(2*((u[i, j+1] - u[i, j-1])*(v[i+1,j] - v[i-1, j])/(4*dx*dy)))-((v[i,j+1]-v[i,j-1])*(v[i,j+1]-v[i, j-1])/(4*dy^2)))


        end
    end

end