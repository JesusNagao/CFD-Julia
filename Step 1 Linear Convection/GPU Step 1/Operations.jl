function kernel_step_1!(u, u_old, c, dt, dx)
    i = threadIdx().x
    #u[i] = u_old[i] - c * dt / dx * (u_old[i] - u_old[i-1])
    @cuprintln(- c * dt / dx)
    return
end