function run_gpu(u::CuArray, v::CuArray, un::CuArray, vn::CuArray, p::CuArray, pn::CuArray, b::CuArray, nx::Int64, ny::Int64, dx::Float64, dy::Float64, rho::Int64, dt::Float64, nit::Int64, F::Int64)

     @cuda blocks=nx threads=ny build_up_b!(b, dx, dy, nx, ny, rho, dt)

     #=for i in range(1, stop=nit)
          un = copy(u)
          vn = copy(v)

          @cuda blocks=nx threads=ny build_up_b!(b, dx, dy, nx, ny, rho, dt)
          pressure_poisson_periodic(p, pn, dx, dy, nx, ny, nit)
     end
     =#
end

function build_up_b!(b, dx::Float64, dy::Float64, nx::Int64, ny::Int64, rho::Int64, dt::Float64)

    i = threadIdx().x
    j = blockIdx().x

    if (i>1 && i<nx && j>1 && j<ny )
     b[i+1, j+1] = (rho * (1 / dt * ((u[i+1, j] - u[i+1, j-1]) / (2 * dx) +
                                        (v[i, j+1] - v[j-1, j+1]) / (2 * dy)) -
                              ((u[i+1, j] - u[i+1, j-1]) / (2 * dx))^2 -
                              2 * ((u[i, j+1] - u[j-1, j+1]) / (2 * dy) *
                                   (v[i+1, j] - v[i+1, j-1]) / (2 * dx))-
                              ((v[i, j+1] - v[j-1, j+1]) / (2 * dy))^2))
     
     # Periodic BC Pressure @ x = 2
     b[i+1, ny] = (rho * (1 / dt * ((u[i+1, 2] - u[i+1, ny-1]) / (2 * dx) +
                                        (v[i, ny] - v[j-1, ny]) / (2 * dy)) -
                              ((u[i+1, 2] - u[i+1, ny-1]) / (2 * dx))^2 -
                              2 * ((u[i, ny] - u[j-1, ny]) / (2 * dy) *
                                   (v[i+1, 2] - v[i+1, ny-1]) / (2 * dx)) -
                              ((v[i, ny] - v[j-1, ny]) / (2 * dy))^2))

     # Periodic BC Pressure @ x = 0
     b[i+1, 2] = (rho * (1 / dt * ((u[i+1, 2] - u[i+1, ny]) / (2 * dx) +
                                        (v[i, 2] - v[j-1, 2]) / (2 * dy)) -
                              ((u[i+1, 2] - u[i+1, ny]) / (2 * dx))^2 -
                              2 * ((u[i, 2] - u[j-1, 2]) / (2 * dy) *
                                   (v[i+1, 2] - v[i+1, ny]) / (2 * dx))-
                              ((v[i, 2] - v[j-1, 2]) / (2 * dy))^2))

    end
    
    return

end

function pressure_poisson_periodic(p::CuArray, pn::CuArray, dx::Float64, dy::Float64, nx::Int64, ny::Int64, nit::Int64)

     for q in range(1, stop=nit)

          @cuda blocks=nx threads=ny calculate_pressure!(p, pn, dx, dy, nx, ny, b)

     end

end

function calculate_pressure!(p::CuArray, pn::CuArray, dx::Float64, dy::Float64, nx::Int64, ny::Int64, b::CuArray)

     i = threadIdx().x
     j = blockIdx().x

     if (i>1 && i<nx && j>1 && j<ny )

          p[i+1, j+1] = (((pn[i+1, j] + pn[i+1, j-1]) * dy^2 +
                              (pn[i, j+1] + pn[i-1, j+1]) * dx^2) /
                              (2 * (dx^2 + dy^2)) -
                              dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[i+1, j+1])

          p[i+1, ny] = (((pn[i+1, 1] + pn[i+1,  ny-1])* dy^2 +
                              (pn[i, ny] + pn[i-1, ny]) * dx^2) /
                         (2 * (dx^2 + dy^2)) -
                         dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[i+1, ny])

          p[i+1, 1] = (((pn[i+1, 2] + pn[i+1, ny])* dy^2 +
                         (pn[i, 1] + pn[i-1, 1]) * dx^2) /
                         (2 * (dx^2 + dy^2)) -
                         dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[i+1, 1])

          
          p[nx, j] = p[nx-1, j]
          p[1, j] = p[2, j]

     end


end