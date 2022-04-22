#using GLMakie

function build_up_b(rho::Int64, dt::Float64, dx::Float64, dy::Float64, u::Matrix{Float64}, v::Matrix{Float64}, b::Matrix{Float64}, nx::Int64, ny::Int64)

    b[2:nx-1,2:ny-1] = (rho * (1 / dt * ((u[2:nx-1, 3:ny] - u[2:nx-1, 1:ny-2]) / (2 * dx) + (v[3:nx, 2:ny-1] - v[1:nx-2, 2:ny-1]) / (2 * dy) - ((u[2:nx-1, 3:ny] - u[2:nx-1, 1:ny-2]) / (2 * dx)) ^ 2 - 2 * ((u[3:nx, 2:ny-1] - u[1:nx-2, 2:ny-1]) / (2*dy) * (v[2:nx-1, 3:ny] - v[2:nx-1, 1:ny-2]) / (2 * dx)) - ((v[3:nx, 2:ny-1] - v[1:nx-2, 2:ny-1]) / (2 * dy))^2)))
    b[2:nx-1,ny-1] = (rho * (1 / dt * ((u[2:nx-1, 1] - u[2:nx-1, ny-2]) / (2 * dx) + (v[3:nx, ny-1] - v[1:nx-2, ny-1]) / (2 * dy) - ((u[2:nx-1, 1] - u[2:nx-1, ny-2]) / (2 * dx)) .^ 2 - 2 * ((u[3:nx, ny-1] - u[1:nx-2, ny-1]) / (2*dy) .* (v[2:nx-1, 1] - v[2:nx-1, ny-2]) / (2 * dx)) - ((v[3:nx, ny-1] - v[1:nx-2, ny-1]) / (2 * dy)).^2)))
    b[2:nx-1, 1] = (rho * (1 / dt * ((u[2:nx-1, 2] - u[2:nx-1, ny-1]) / (2 * dx) + (v[3:nx, 1] - v[1:nx-2, 1]) / (2 * dy) - ((u[2:nx-1, 2] - u[2:nx-1, ny-1]) / (2 * dx)) .^ 2 - 2 * ((u[3:nx, 1] - u[1:nx-2, ny-1]) / (2*dy) .* (v[2:nx-1, 1] - v[2:nx-1, ny-1]) / (2 * dx)) - ((v[3:nx, 1] - v[1:nx-2, 1]) / (2 * dy)).^2)))
    
end

function pressure_poisson_periodic(p::Matrix{Float64}, dx::Float64, dy::Float64, nit::Int64, nx::Int64, ny::Int64)

    for j in range(1, stop=nit)
    
        pn = Matrix{Float64}(undef, nx, ny)
        pn = p

        p[2:nx-1, 2:ny-1] = (((pn[2:nx-1, 3:ny] + pn[2:nx-1, 1:ny-2]) * dy^2 + (pn[3:nx, 2:ny-1] + pn[1:nx-2, 2:ny-1]) * dx^2) / (2 * (dx^2 + dy^2)) - dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:nx-1, 2:nx-1])
        p[2:nx-1, ny-1] = (((pn[2:nx-1, 1] + pn[2:nx-1, ny-2])* dy^2 + (pn[3:nx, ny-1] + pn[1:nx-2, ny-1]) * dx^2) / (2 * (dx^2 + dy^2)) - dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:nx-1, ny-1])
        p[2:nx-1, 1] = (((pn[2:nx-1, 2] + pn[2:nx-1, ny-1])* dy^2 + (pn[3:nx, 1] + pn[1:ny-2, 1]) * dx^2) / (2 * (dx^2 + dy^2)) - dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:nx-1, 1])
        p[nx-1, :] =p[nx-2, :]
        p[1, :] = p[2, :]

    end

end