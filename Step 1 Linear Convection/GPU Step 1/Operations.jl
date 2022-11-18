function kernel_step_1!(u_new, uf, ub, c, dt, dx)
    i = threadIdx().x

    u_new[i] = uf[i] - c * dt / dx * (uf[i] - ub[i])

    #@cuprint(u_new)
    return
end