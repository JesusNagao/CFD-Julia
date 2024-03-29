function build_up_b(u, v, b, nx, ny, dx, dy, rho)

    
    b[2:nx-1, 2:nx-1] = (rho .* (1 / dt .* ((u[2:nx-1, 3:nx] - u[2:nx-1, 1:nx-2]) / (2 .* dx) +
                                      (v[3:nx, 2:nx-1] - v[1:nx-2, 2:nx-1]) / (2 .* dy)) -
                            ((u[2:nx-1, 3:nx] - u[2:nx-1, 1:nx-2]) / (2 .* dx)) .^ 2 -
                            2 .* ((u[3:nx, 2:nx-1] - u[1:nx-2, 2:nx-1]) / (2 .* dy) .*
                                 (v[2:nx-1, 3:nx] - v[2:nx-1, 1:nx-2]) / (2 .* dx))-
                            ((v[3:nx, 2:nx-1] - v[1:nx-2, 2:nx-1]) / (2 .* dy)) .^ 2))
    
    # Periodic BC Pressure @ x = 2
    b[2:nx-1, nx-1] = (rho .* (1 / dt .* ((u[2:nx-1, 1] - u[2:nx-1,nx-2]) / (2 .* dx) +
                                    (v[3:nx, nx-1] - v[1:nx-2, nx-1]) / (2 .* dy)) -
                          ((u[2:nx-1, 1] - u[2:nx-1, nx-2]) / (2 .* dx)) .^ 2 -
                          2 .* ((u[3:nx, nx-1] - u[1:nx-2, nx-1]) / (2 .* dy) .*
                               (v[2:nx-1, 1] - v[2:nx-1, nx-2]) / (2 .* dx)) -
                          ((v[3:nx, nx-1] - v[1:nx-2, nx-1]) / (2 .* dy)) .^ 2))

    # Periodic BC Pressure @ x = 0
    b[2:nx-1, 1] = (rho .* (1 / dt .* ((u[2:nx-1, 1] - u[2:nx-1, nx-1]) / (2 .* dx) +
                                   (v[3:nx, 1] - v[1:nx-2, 1]) / (2 .* dy)) -
                         ((u[2:nx-1, 1] - u[2:nx-1, nx-1]) / (2 .* dx)) .^ 2 -
                         2 .* ((u[3:nx, 1] - u[1:nx-2, 1]) / (2 .* dy) .*
                              (v[2:nx-1, 1] - v[2:nx-1, nx-1]) / (2 .* dx))-
                         ((v[3:nx, 1] - v[1:nx-2, 1]) / (2 .* dy)) .^ 2))

end

function pressure_poisson_periodic(pn, p, nx, ny, dx, dy, b, nit)
    
    for q in range(1, stop=nit)
        pn[:,:] = p[:,:]
        p[2:nx-1, 2:nx-1] = (((pn[2:nx-1, 3:nx] + pn[2:nx-1, 1:nx-2]) * dy^2 +
                          (pn[3:nx, 2:nx-1] + pn[1:nx-2, 2:nx-1]) * dx^2) /
                         (2 * (dx^2 + dy^2)) -
                         dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:nx-1, 2:nx-1])
    
        # Periodic BC Pressure @ x = 2
        p[2:nx-1, nx-1] = (((pn[2:nx-1, 1] + pn[2:nx-1, nx-2])* dy^2 +
                        (pn[3:nx, nx-1] + pn[1:nx-2, nx-1]) * dx^2) /
                       (2 * (dx^2 + dy^2)) -
                       dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:nx-1, nx-1])
    
        # Periodic BC Pressure @ x = 0
        p[2:nx-1, 1] = (((pn[2:nx-1, 2] + pn[2:nx-1, nx-1])* dy^2 +
                       (pn[3:nx, 1] + pn[1:nx-2, 1]) * dx^2) /
                      (2 * (dx^2 + dy^2)) -
                      dx^2 * dy^2 / (2 * (dx^2 + dy^2)) * b[2:nx-1, 1])
        
        # Wall boundary conditions, pressure
        p[nx-1, :] = p[nx-2, :]  # dp/dy = 0 at y = 2
        p[1, :] = p[2, :]  # dp/dy = 0 at y = 0
    end

end

function run(u, v, un, vn, nx, ny, dx, dy, rho, F, b, nit, pn, p, nu)

    f = Figure(resolution = (800, 800))
    Axis(f[1,1], backgroundcolor = "black")
    
    for i in range(0, stop=1)
    #record(f, "Channel Flow.gif", 1:499) do i
        un[:,:] = u[:,:]
        vn[:,:] = v[:,:]
    
        build_up_b(u, v, b, nx, ny, dx, dy, rho)
        pressure_poisson_periodic(pn, p, nx, ny, dx, dy, b, nit)
    
        u[2:nx-1, 2:nx-1] = (un[2:nx-1, 2:nx-1] -
                         un[2:nx-1, 2:nx-1] .* dt / dx .* 
                        (un[2:nx-1, 2:nx-1] - un[2:nx-1, 1:nx-2]) -
                         vn[2:nx-1, 2:nx-1] .* dt / dy .* 
                        (un[2:nx-1, 2:nx-1] - un[1:nx-2, 2:nx-1]) -
                         dt / (2 .* rho .* dx) .* 
                        (p[2:nx-1, 3:nx] - p[2:nx-1, 1:nx-2]) .+
                         nu .* (dt / dx^2 .* 
                        (un[2:nx-1, 3:nx] - 2 .* un[2:nx-1, 2:nx-1] .+ un[2:nx-1, 1:nx-2]) .+
                         dt / dy^2 .* 
                        (un[3:nx, 2:nx-1] - 2 .* un[2:nx-1, 2:nx-1] .+ un[1:nx-2, 2:nx-1])) .+ 
                         F .* dt)
    
        v[2:nx-1, 2:nx-1] = (vn[2:nx-1, 2:nx-1] -
                         un[2:nx-1, 2:nx-1] .* dt / dx .* 
                        (vn[2:nx-1, 2:nx-1] - vn[2:nx-1, 1:nx-2]) -
                         vn[2:nx-1, 2:nx-1] .* dt / dy .* 
                        (vn[2:nx-1, 2:nx-1] - vn[1:nx-2, 2:nx-1]) -
                         dt / (2 .* rho .* dy) .* 
                        (p[3:nx, 2:nx-1] - p[1:nx-2, 2:nx-1]) .+
                         nu .* (dt / dx^2 .*
                        (vn[2:nx-1, 3:nx] - 2 .* vn[2:nx-1, 2:nx-1] .+ vn[2:nx-1, 1:nx-2]) .+
                         dt / dy^2 .* 
                        (vn[3:nx, 2:nx-1] - 2 .* vn[2:nx-1, 2:nx-1] .+ vn[1:nx-2, 2:nx-1])))
    
        # Periodic BC u @ x = 2     
        u[2:nx-1, nx-1] = (un[2:nx-1, nx-1] - un[2:nx-1, nx-1] .* dt / dx .* 
                      (un[2:nx-1, nx-1] - un[2:nx-1, nx-2]) -
                       vn[2:nx-1, nx-1] .* dt / dy .* 
                      (un[2:nx-1, nx-1] - un[1:nx-2, nx-1]) -
                       dt / (2 .* rho .* dx) .*
                      (p[2:nx-1, 1] - p[2:nx-1, nx-2]) .+ 
                       nu .* (dt / dx^2 .* 
                      (un[2:nx-1, 1] - 2 .* un[2:nx-1, nx-1] .+ un[2:nx-1, nx-2]) .+
                       dt / dy^2 .* 
                      (un[3:nx, nx-1] - 2 .* un[2:nx-1, nx-1] .+ un[1:nx-2, nx-1])) .+ F .* dt)
    
        # Periodic BC u @ x = 0
        u[2:nx-1, 1] = (un[2:nx-1, 1] - un[2:nx-1, 1] .* dt / dx .*
                     (un[2:nx-1, 1] - un[2:nx-1, nx-1]) -
                      vn[2:nx-1, 1] .* dt / dy .* 
                     (un[2:nx-1, 1] - un[1:nx-2, 1]) - 
                      dt / (2 .* rho .* dx) .* 
                     (p[2:nx-1, 2] - p[2:nx-1, nx-1]) .+ 
                      nu .* (dt / dx^2 .* 
                     (un[2:nx-1, 2] - 2 .* un[2:nx-1, 1] .+ un[2:nx-1, nx-1]) .+
                      dt / dy^2 .*
                     (un[3:nx, 1] - 2 .* un[2:nx-1, 1] .+ un[1:nx-2, 1])) .+ F .* dt)
    
        # Periodic BC v @ x = 2
        v[2:nx-1, nx-1] = (vn[2:nx-1, nx-1] - un[2:nx-1, nx-1] .* dt / dx .*
                      (vn[2:nx-1, nx-1] - vn[2:nx-1, nx-2]) - 
                       vn[2:nx-1, nx-1] .* dt / dy .*
                      (vn[2:nx-1, nx-1] - vn[1:nx-2, nx-1]) -
                       dt / (2 .* rho .* dy) .* 
                      (p[3:nx, nx-1] - p[1:nx-2, nx-1]) .+
                       nu .* (dt / dx^2 .*
                      (vn[2:nx-1, 1] - 2 .* vn[2:nx-1, nx-1] .+ vn[2:nx-1, nx-2]) .+
                       dt / dy^2 .*
                      (vn[3:nx, nx-1] - 2 .* vn[2:nx-1, nx-1] .+ vn[1:nx-2, nx-1])))
    
        # Periodic BC v @ x = 0
        v[2:nx-1, 1] = (vn[2:nx-1, 1] - un[2:nx-1, 1] .* dt / dx .*
                     (vn[2:nx-1, 1] - vn[2:nx-1, nx-1]) -
                      vn[2:nx-1, 1] .* dt / dy .*
                     (vn[2:nx-1, 1] - vn[1:nx-2, 1]) -
                      dt / (2 .* rho .* dy) .* 
                     (p[3:nx, 1] - p[1:nx-2, 1]) .+
                      nu .* (dt / dx^2 .* 
                     (vn[2:nx-1, 2] - 2 .* vn[2:nx-1, 1] .+ vn[2:nx-1, nx-1]) .+
                      dt / dy^2 .* 
                     (vn[3:nx, 1] - 2 .* vn[2:nx-1, 1] .+ vn[1:nx-2, 1])))
    
    
        # Wall BC: u,v = 0 @ y = 0,2
        u[1, :] .= 0.0
        u[nx-1, :] .= 0.0
        v[1, :] .= 0.0
        v[nx-1, :] .= 0.0
        
        strength = vec(sqrt.(u[2:nx-2, 2:nx-2] .^ 2 .+ v[2:nx-2, 2:nx-2] .^ 2))
        #println(length(strength))
        arrows!(x[2:nx-2], y[2:nx-2], v[2:nx-2, 2:nx-2], u[2:nx-2, 2:nx-2], lengthscale = 0.01, arrowcolor = strength, linecolor = strength)
    end
end

function run_kernel()
    



end