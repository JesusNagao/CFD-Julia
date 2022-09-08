using GLMakie

const nx = 41
const ny = 41
const nz = 41
const nt = 10
const nit = 51
const c = 1
const dx = 2 / (nx - 1)
const dy = 2 / (ny - 1)
const dz = 2 / (nz - 1)

x = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = Array{Float64}([i for i in range(0.0, stop=2.0, step=dy)])
z = Array{Float64}([i for i in range(0.0, stop=2.0, step=dz)])

const rho = 1
const nu = 0.1
const dt = 0.01

u = zeros(nx, ny, nz)
v = zeros(nx, ny, nz)
w = zeros(nx, ny, nz)
p = zeros(nx, ny, nz)
un = zeros(nx, ny, nz)
vn = zeros(nx, ny, nz)
wn = zeros(nx, ny, nz)
pn = zeros(nx, ny, nz)

function run(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, p::Array{Float64}, un::Array{Float64}, vn::Array{Float64}, wn::Array{Float64}, pn::Array{Float64}, nit::Int64, nx::Int64, ny::Int64, nz::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64)

    #f = Figure(resolution = (800, 800))
    #Axis(f[1,1], backgroundcolor = "black")
    n=0

    while (n<nit)
    #record(f, "Channel Flow.gif", 1:499) do i
        for i in range(2, stop=nx-2)
            for j in range(2, stop=ny-2)
                for k in range(2, stop=nz-2)
                    un[i,j,k] = calc_u(u, i, j, k, rho, dx, dy, dz, dt)
                    vn[i,j,k] = calc_v(v, i, j, k, rho, dx, dy, dz, dt)
                    wn[i,j,k] = calc_w(w, i, j, k, rho, dx, dy, dz, dt)
                    pn[i,j,k] = calc_p(u, v, w, dx, dy, dz, rho, i, j, k)
                end
            end
        end

        u[:,:,:] = un[:,:,:]
        v[:,:,:] = vn[:,:,:]
        w[:,:,:] = wn[:,:,:]
        p[:,:,:] = pn[:,:,:]
        

        u[1, :, :] .= 0.0
        u[nx-1, :, :] .= 0.0
        v[1, :, :] .= 0.0
        v[ny-1, :, :] .= 0.0
        w[1, :, :] .= 0.0
        w[nz-1, :, :] .= 0.0

        u[:, 1, :] .= 0.0
        u[:, nx-1, :] .= 0.0
        v[:, 1, :] .= 0.0
        v[:, ny-1, :] .= 0.0
        w[:, 1, :] .= 0.0
        w[:, nz-1, :] .= 0.0

        
        periodic_bc(u, v, w, p, nx, ny, nz, dx, dy, dz, dt, rho)
        
        #strength = vec(sqrt.(u[2:nx-2, 2:ny-2, 2:nz-2] .^ 2 .+ v[2:nx-2, 2:ny-2, 2:nz-2] .^ 2 .+ w[2:nx-2, 2:ny-2, 2:nz-2] .^ 2))
        #arrows!(x[2:nx-2], y[2:ny-2], z[2:nz-2], w[2:nx-2, 2:ny-2, 2:nz-2], v[2:nx-2, 2:ny-2, 2:nz-2], u[2:nx-2, 2:ny-2, 2:nz-2], lengthscale = 0.01, arrowcolor = strength, linecolor = strength)

        n = n + 1

    end

end

function calc_u(u::Array{Float64}, i::Int64, j::Int64, k::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64)

    un = u[i,j,k]-(p[i+1,j,k]-p[i-1,j,k])/(2*dx*rho)
    + nu*(((u[i+1,j,k]-2*u[i,j,k]+u[i-1,j,k])/(dx^2))+((u[i,j+1,k]-2*u[i,j,k]+u[i,j-1,k])/(dy^2))+((u[i,j,k+1]-2*u[i,j,k]+u[i,j,k-1])/(dz^2)))
    - u[i,j,k]*((u[i,j,k]-u[i-1,j,k])/(dx)) - v[i,j,k]*((u[i,j,k]-u[i,j-1,k])/(dy)) - w[i,j,k]*((u[i,j,k]-u[i,j,k-1])/(dz))*dt

    return un
end

function calc_v(v::Array{Float64}, i::Int64, j::Int64, k::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64)

    vn = v[i,j,k]-(p[i,j+1,k]-p[i,j-1,k])/(2*dy*rho)
    + nu*(((v[i+1,j,k]-2*v[i,j,k]+v[i-1,j,k])/(dx^2))+(v[i,j+1,k]-2*v[i,j,k]+v[i,j-1,k])/(dy^2))+((v[i,j,k+1]-2*v[i,j,k]+v[i,j,k-1])/(dz^2))
    - u[i,j,k]*((v[i,j,k]-v[i-1,j,k])/(dx)) - v[i,j,k]*((v[i,j,k]-v[i,j-1,k])/(dy)) - w[i,j,k]*((v[i,j,k]-v[i,j,k-1])/(dz))*dt

    return vn
end

function calc_w(w::Array{Float64}, i::Int64, j::Int64, k::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64)

    wn = w[i,j,k]-(p[i,j+1,k]-p[i,j-1,k])/(2*dz*rho)
    + nu*(((w[i+1,j,k]-2*w[i,j,k]+w[i-1,j,k])/(dx^2))+(w[i,j+1,k]-2*w[i,j,k]+w[i,j-1,k])/(dy^2))+((w[i,j,k+1]-2*w[i,j,k]+w[i,j,k-1])/(dz^2))
    - u[i,j,k]*((w[i,j,k]-w[i-1,j,k])/(dx)) - v[i,j,k]*((w[i,j,k]-w[i,j-1,k])/(dy)) - w[i,j,k]*((w[i,j,k]-w[i,j,k-1])/(dz))*dt

    return wn
end

function calc_p(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, dx::Float64, dy::Float64, dz::Float64, rho::Int64, i::Int64, j::Int64, k::Int64)

    p = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) *
    (((u[i+1,j,k]-u[i-1,j,k])/(2*dx))^2 + ((v[i,j-1,k]-v[i,j+1,k])/(2*dy))^2 + ((w[i,j,k-1]-w[i,j,k-1])/(2*dz))^2
    + 2*((u[i,j+1,k]-u[i,j-1,k])/(2*dy))*((v[i+1,j,k]-v[i-1,j,k])/(2*dx)) + 2*((w[i,j+1,k]-w[i,j-1,k])/(2*dy))*((v[i,j,k+1]-v[i,j,k-1])/(2*dz))
    + 2*((u[i,j,k+1]-u[i,j,k-1])/(2*dz))*((w[i+1,j,k]-w[i-1,j,k])/(2*dx)))

    return p
end

function periodic_bc(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, p::Array{Float64}, nx::Int64, ny::Int64, nz::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64, rho::Int64)

    #Periodic BC at z=0 for u
    u[2:nx-1,2:ny-1,1] = u[2:nx-1,2:ny-1,1]-(p[3:nx,2:ny-1,1]-p[1:nx-2,2:ny-1,1])/(2*dx*rho)
    + nu*(((u[3:nx,2:ny-1,1]-2*u[2:nx-1,2:ny-1,1]+u[1:nx-2,2:ny-1,1])/(dx^2))+((u[2:nx-1,3:ny,1]-2*u[2:nx-1,3:ny,1]+u[2:nx-1,1:ny-2,1])/(dy^2))+((u[2:nx-1,2:ny-1,2]-2*u[2:nx-1,2:ny-1,1]+u[2:nx-1,2:ny-1, nz-1])/(dz^2)))
    - u[2:nx-1,2:ny-1,1]*((u[2:nx-1,2:ny-1,1]-u[1:nx-2,2:ny-1,1])/(dx)) - v[2:nx-1,2:ny-1,1]*((u[2:nx-1,2:ny-1,1]-u[2:nx-1,1:ny-2,1])/(dy)) - w[2:nx-1,2:ny-1,1]*((u[2:nx-1,2:ny-1,1]-u[2:nx-1,2:ny-1,nz-1])/(dz))*dt

    #Periodic BC at z=2 for u
    u[2:nx-1,2:ny-1,nz-1] = u[2:nx-1,2:ny-1,nz-1]-(p[3:nx,2:ny-1,nz-1]-p[1:nx-2,2:ny-1,nz-1])/(2*dx*rho)
    + nu*(((u[3:nx,2:ny-1,nz-1]-2*u[2:nx-1,2:ny-1,nz-1]+u[1:nx-2,2:ny-1,nz-1])/(dx^2))+((u[2:nx-1,3:ny,nz-1]-2*u[2:nx-1,3:ny,nz-1]+u[2:nx-1,1:ny-2,nz-1])/(dy^2))+((u[2:nx-1,2:ny-1, 1]-2*u[2:nx-1,2:ny-1,nz-1]+u[2:nx-1,2:ny-1,nz-2])/(dz^2)))
    - u[2:nx-1,2:ny-1,nz-1]*((u[2:nx-1,2:ny-1,nz-1]-u[1:nx-2,2:ny-1,nz-1])/(dx)) - v[2:nx-1,2:ny-1,nz-1]*((u[2:nx-1,2:ny-1,nz-1]-u[2:nx-1,1:ny-2,nz-1])/(dy)) - w[2:nx-1,2:ny-1,nz-1]*((u[2:nx-1,2:ny-1,nz-1]-u[2:nx-1,2:ny-1,nz-2])/(dz))*dt
    
    #Periodic BC at z=0 for v
    v[2:nx-1,2:ny-1,1] = v[2:nx-1,2:ny-1,1]-(p[2:nx-1,3:ny,1]-p[2:nx-1,1:ny-2,1])/(2*dy*rho)
    + nu*(((v[3:nx,2:ny-1,1]-2*v[2:nx-1,2:ny-1,1]+v[1:nx-2,2:ny-1,1])/(dx^2))+(v[2:nx-1,3:ny,1]-2*v[2:nx-1,2:ny-1,1]+v[2:nx-1,1:ny-2,1])/(dy^2))+((v[2:nx-1,2:ny-1,2]-2*v[2:nx-1,2:ny-1,1]+v[2:nx-1,2:ny-1,nz-1])/(dz^2))
    - u[2:nx-1,2:ny-1,1]*((v[2:nx-1,2:ny-1,1]-v[1:nx-2,2:ny-1,1])/(dx)) - v[2:nx-1,2:ny-1,1]*((v[2:nx-1,2:ny-1,1]-v[2:nx-1,1:ny-2,1])/(dy)) - w[2:nx-1,2:ny-1,1]*((v[2:nx-1,2:ny-1,1]-v[2:nx-1,2:ny-1,nx-1])/(dz))*dt

    #Periodic BC at z=2 for v
    v[2:nx-1,2:ny-1,nz-1] = v[2:nx-1,2:ny-1,nz-1]-(p[2:nx-1,3:ny,nz-1]-p[2:nx-1,1:ny-2,nz-1])/(2*dy*rho)
    + nu*(((v[3:nx,2:ny-1,nz-1]-2*v[2:nx-1,2:ny-1,nz-1]+v[1:nx-2,2:ny-1,nz-1])/(dx^2))+(v[2:nx-1,3:ny,nz-1]-2*v[2:nx-1,2:ny-1,nz-1]+v[2:nx-1,1:ny-2,nz-1])/(dy^2))+((v[2:nx-1,2:ny-1,1]-2*v[2:nx-1,2:ny-1,nz-1]+v[2:nx-1,2:ny-1,nz-2])/(dz^2))
    - u[2:nx-1,2:ny-1,nz-1]*((v[2:nx-1,2:ny-1,nz-1]-v[1:nx-2,2:ny-1,nz-1])/(dx)) - v[2:nx-1,2:ny-1,nz-1]*((v[2:nx-1,2:ny-1,nz-1]-v[2:nx-1,1:ny-2,nz-1])/(dy)) - w[2:nx-1,2:ny-1,nz-1]*((v[2:nx-1,2:ny-1,nz-1]-v[2:nx-1,2:ny-1,nz-2])/(dz))*dt
    
    #Periodic BC at z=0 for w
    w[2:nx-1,2:ny-1,1] = w[2:nx-1,2:ny-1,1]-(p[2:nx-1,3:ny,1]-p[2:nx-1,1:ny-2,1])/(2*dz*rho)
    + nu*(((w[3:nx,2:ny-1,1]-2*w[2:nx-1,2:ny-1,1]+w[1:nx-2,2:ny-1,1])/(dx^2))+(w[2:nx-1,3:ny,1]-2*w[2:nx-1,2:ny-1,1]+w[2:nx-1,1:ny-2,1])/(dy^2))+((w[2:nx-1,2:ny-1,2]-2*w[2:nx-1,2:ny-1,1]+w[2:nx-1,2:ny-1,nz-1])/(dz^2))
    - u[2:nx-1,2:ny-1,1]*((w[2:nx-1,2:ny-1,1]-w[1:nx-2,2:ny-1,1])/(dx)) - v[2:nx-1,2:ny-1,1]*((w[2:nx-1,2:ny-1,1]-w[2:nx-1,1:ny-2,1])/(dy)) - w[2:nx-1,2:ny-1,1]*((w[2:nx-1,2:ny-1,1]-w[2:nx-1,2:ny-1,nz-1])/(dz))*dt

    #Periodic BC at z=2 for w
    w[2:nx-1,2:ny-1,nz-1] = w[2:nx-1,2:ny-1,nz-1]-(p[2:nx-1,3:ny,nz-1]-p[2:nx-1,1:ny-2,nz-1])/(2*dz*rho)
    + nu*(((w[3:nx,2:ny-1,nz-1]-2*w[2:nx-1,2:ny-1,nz-1]+w[1:nx-2,2:ny-1,nz-1])/(dx^2))+(w[2:nx-1,3:ny,nz-1]-2*w[2:nx-1,2:ny-1,nz-1]+w[2:nx-1,1:ny-2,nz-1])/(dy^2))+((w[2:nx-1,2:ny-1,1]-2*w[2:nx-1,2:ny-1,nz-1]+w[2:nx-1,2:ny-1,nz-2])/(dz^2))
    - u[2:nx-1,2:ny-1,nz-1]*((w[2:nx-1,2:ny-1,nz-1]-w[1:nx-2,2:ny-1,nz-1])/(dx)) - v[2:nx-1,2:ny-1,nz-1]*((w[2:nx-1,2:ny-1,nz-1]-w[2:nx-1,1:ny-2,nz-1])/(dy)) - w[2:nx-1,2:ny-1,nz-1]*((w[2:nx-1,2:ny-1,nz-1]-w[2:nx-1,2:ny-1,nz-2])/(dz))*dt

    #Periodic BC at z=0 for p
    p[2:nx-1,2:ny-1,1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) *
    (((u[3:nx,2:ny-1,1]-u[1:nx-2,2:ny-1,1])/(2*dx))^2 + ((v[2:nx-1,1:ny-2,1]-v[2:nx-1,3:ny,1])/(2*dy))^2 + ((w[2:nx-1,2:ny-1,nz-1]-w[2:nx-1,2:ny-1,nz-1])/(2*dz))^2
    + 2*((u[2:nx-1,3:ny,1]-u[2:nx-1,1:ny-2,1])/(2*dy))*((v[3:nx,2:ny-1,1]-v[1:nx-2,2:ny-1,1])/(2*dx)) + 2*((w[2:nx-1,3:ny,1]-w[2:nx-1,1:ny-2,1])/(2*dy))*((v[2:nx-1,2:ny-1,2]-v[2:nx-1,2:ny-1,nz-1])/(2*dz))
    + 2*((u[2:nx-1,2:ny-1,2]-u[2:nx-1,2:ny-1,nz-1])/(2*dz))*((w[3:nx,2:ny-1,1]-w[1:nx-2,2:ny-1,1])/(2*dx)))
    
    #Periodic BC at z=2 for p
    p[2:nx-1,2:ny-1,nz-1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) *
    (((u[3:nx,2:ny-1,nz-1]-u[1:nx-2,2:ny-1,nz-1])/(2*dx))^2 + ((v[2:nx-1,1:ny-2,nz-1]-v[2:nx-1,3:ny,nz-1])/(2*dy))^2 + ((w[2:nx-1,2:ny-1,nz-2]-w[2:nx-1,2:ny-1,nz-2])/(2*dz))^2
    + 2*((u[2:nx-1,3:ny,nz-1]-u[2:nx-1,1:ny-2,nz-1])/(2*dy))*((v[3:nx,2:ny-1,nz-1]-v[1:nx-2,2:ny-1,nz-1])/(2*dx)) + 2*((w[2:nx-1,3:ny,nz-1]-w[2:nx-1,1:ny-2,nz-1])/(2*dy))*((v[2:nx-1,2:ny-1,1]-v[2:nx-1,2:ny-1,nz-2])/(2*dz))
    + 2*((u[2:nx-1,2:ny-1,1]-u[2:nx-1,2:ny-1,nz-2])/(2*dz))*((w[3:nx,2:ny-1,nz-1]-w[1:nx-2,2:ny-1,nz-1])/(2*dx)))


    p[nx-1, ny-1,:] = p[nx-2, ny-2,:]
    p[1, 1, :] = p[2, 2, :] 

end

run(u, v, w, p, un, vn, wn, pn, 200, nx, ny, nz, rho, dx, dy, dz, dt)