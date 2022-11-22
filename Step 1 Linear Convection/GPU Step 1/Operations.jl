function kernel_step_1!(u, u_old, c, dt, dx, nx)
    i = threadIdx().x

    if (i>1 && i<nx)
        u[i] = u_old[i] - c * dt / dx * (u_old[i] - u_old[i-1])
    end

    return
end

function initialize!(u, nx)
    i = threadIdx().x

    if (i>Int(round(nx/3)) && i<Int(round(2*nx/3)))

        u[i] = 2.0

    end

    return

end