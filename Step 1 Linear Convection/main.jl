using Plots

const nx = 41
const dx = 2 / (nx-1)
const nt = 25
const dt = 0.025
const c = 1

u = Array{Float64}(undef , nx);
u = ones(nx)

for i in range(Int(round(nx/3)), stop=Int(round(2*nx/3)))

    u[i] = 2;

end

pb = plot(u)

for i in range(1, stop=nt)
    u_old = Array{Float64}(undef, nx)
    u_old = copy(u)
    for j in range(2, stop=nx)

        u[j] = u_old[j] - c * dt / dx * (u_old[j]- u_old[j-1])

    end
    
end

pa = plot(u)

plot(pb, pa, layout=(2,1))
