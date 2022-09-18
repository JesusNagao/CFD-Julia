using GLMakie

function run(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, p::Array{Float64}, un::Array{Float64}, vn::Array{Float64}, wn::Array{Float64}, pn::Array{Float64}, nx::Int64, ny::Int64, nz::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64, F::Int64, X::Array{Float64}, Y::Array{Float64}, Z::Array{Float64})
    

    u[:,:,:] = un[:,:,:]
    v[:,:,:] = vn[:,:,:]
    w[:,:,:] = wn[:,:,:]
    p[:,:,:] = pn[:,:,:]

    pn[nx-1, :,:] = pn[nx-2, :,:]
    pn[1, :, :] = pn[2, :, :]
    pn[:, ny-1,:] = pn[:, ny-2,:]
    pn[:, 1, :] = pn[:, 2, :]

    pn[2:nx-1, 2:ny-1, 2:nz-1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) * (((u[3:nx,2:ny-1,2:nz-1]-u[1:nx-2,2:ny-1,2:nz-1])/(2*dx)).^2 .+ ((v[2:nx-1,1:ny-2,2:nz-1]-v[2:nx-1,3:ny,2:nz-1])/(2*dy)).^2 .+ ((w[2:nx-1,2:ny-1,1:nz-2]-w[2:nx-1,2:ny-1,1:nz-2])/(2*dz)).^2 .+ 2*((u[2:nx-1,3:ny,2:nz-1]-u[2:nx-1,1:ny-2,2:nz-1])/(2*dy)).*((v[3:nx,2:ny-1,2:nz-1]-v[1:nx-2,2:ny-1,2:nz-1])/(2*dx)) + 2*((w[2:nx-1,3:ny,2:nz-1]-w[2:nx-1,1:ny-2,2:nz-1])/(2*dy)).*((v[2:nx-1,2:ny-1,3:nz]-v[2:nx-1,2:ny-1,1:nz-2])/(2*dz)) + 2*((u[2:nx-1,2:ny-1,3:nz]-u[2:nx-1,2:ny-1,1:nz-2])/(2*dz)).*((w[3:nx,2:ny-1,2:nz-1]-w[1:nx-2,2:ny-1,2:nz-1])/(2*dx)))
    un[2:nx-1, 2:ny-1, 2:nz-1] = u[2:nx-1, 2:ny-1, 2:nz-1] .+ (-(p[3:nx, 2:ny-1, 2:nz-1] - p[1:nx-2, 2:ny-1, 2:nz-1])/(2*rho*dx) + nu*((u[3:nx, 2:ny-1, 2:nz-1] - 2*u[2:nx-1, 2:ny-1, 2:nz-1] + u[1:nx-2, 2:ny-1, 2:nz-1])/(dx^2) + (u[2:nx-1, 3:ny, 2:nz-1] - 2*u[2:nx-1, 2:ny-1, 2:nz-1] + u[2:nx-1, 1:ny-2, 2:nz-1])/(dy^2) + (u[2:nx-1, 2:ny-1, 3:nx] - 2*u[2:nx-1, 2:ny-1, 2:nz-1] + u[2:nx-1, 2:ny-1, 1:nz-2])/(dz^2)) - u[2:nx-1, 2:ny-1, 2:nz-1].*(u[2:nx-1, 2:ny-1, 2:nz-1] - u[1:nx-2, 2:ny-1, 2:nz-1])/(dx) - v[2:nx-1, 2:ny-1, 2:nz-1].*(u[2:nx-1, 2:ny-1, 2:nz-1] - u[2:nx-1, 1:ny-2, 2:nz-1])/(dy) - w[2:nx-1, 2:ny-1, 2:nz-1].*(u[2:nx-1, 2:ny-1, 2:nz-1] - u[2:nx-1, 2:ny-1, 1:nz-2])/(dz))*dt .+ F*dt
    vn[2:nx-1, 2:ny-1, 2:nz-1] = v[2:nx-1, 2:ny-1, 2:nz-1] .+ (-(p[2:nx-1, 3:ny, 2:nz-1] - p[2:nx-1, 1:ny-2, 2:nx-1])/(2*dy*rho) + nu*((v[3:nx, 2:ny-1, 2:nz-1] - 2*v[2:nx-1, 2:ny-1, 2:nz-1] + v[1:nx-2, 2:ny-1, 2:nz-1])/(dx^2) + (v[2:nx-1, 3:ny, 2:nz-1] - 2*v[2:nx-1, 2:ny-1, 2:nz-1] + v[2:nx-1, 1:ny-2, 2:nz-1])/(dy^2) + (v[2:nx-1, 2:ny-1, 3:nz] - 2*v[2:nx-1, 2:ny-1, 2:nz-1] + v[2:nx-1, 2:ny-1, 1:nz-2])/(dz^2)) - u[2:nx-1, 2:ny-1, 2:nz-1].*(v[2:nx-1, 2:ny-1, 2:nz-1] - v[1:nx-2, 2:ny-1, 2:nz-1])/(dx) - v[2:nx-1, 2:ny-1, 2:nz-1].*(v[2:nx-1, 2:ny-1, 2:nz-1] - v[2:nx-1, 1:ny-2, 2:nz-1])/(dy) - w[2:nx-1, 2:ny-1, 2:nz-1].*(v[2:nx-1, 2:ny-1, 2:nz-1] - v[2:nx-1, 2:ny-1, 1:nz-2])/(dz))*dt
    wn[2:nx-1, 2:ny-1, 2:nz-1] = w[2:nx-1, 2:ny-1, 2:nz-1] .+ (-(p[2:nx-1, 2:ny-1, 3:nz] - p[2:nx-1, 2:ny-1, 1:nz-2])/(2*dz*rho) + nu*((w[3:nx, 2:ny-1, 2:nz-1] - 2*w[2:nx-1, 2:ny-1, 2:nz-1] + w[1:nx-2, 2:ny-1, 2:nz-1])/(dx^2) + (w[2:nx-1, 3:ny, 2:nz-1] - 2*w[2:nx-1, 2:ny-1, 2:nz-1] + w[2:nx-1, 1:ny-2, 2:nz-1])/(dy^2) + (w[2:nx-1, 2:ny-1, 3:nz] - 2*w[2:nx-1, 2:ny-1, 2:nz-1] + w[2:nx-1, 2:ny-1, 1:nz-2])/(dz^2)) - u[2:nx-1, 2:ny-1, 2:nz-1].*(w[2:nx-1, 2:ny-1, 2:nz-1] - w[1:nx-2, 2:ny-1, 2:nz-1])/(dx) - v[2:nx-1, 2:ny-1, 2:nz-1].*(w[2:nx-1, 2:ny-1, 2:nz-1] - w[2:nx-1, 1:ny-2, 2:nz-1])/(dy) - w[2:nx-1, 2:ny-1, 2:nz-1].*(w[2:nx-1, 2:ny-1, 2:nz-1] - w[2:nx-1, 2:ny-1, 1:nz-2])/(dz))*dt

    un[1, :, :] .= 0.0
    un[nx-1, :, :] .= 0.0
    vn[1, :, :] .= 0.0
    vn[nx-1, :, :] .= 0.0
    wn[1, :, :] .= 0.0
    wn[nx-1, :, :] .= 0.0

    un[:, 1, :] .= 0.0
    un[:, ny-1, :] .= 0.0
    vn[:, 1, :] .= 0.0
    vn[:, ny-1, :] .= 0.0
    wn[:, 1, :] .= 0.0
    wn[:, ny-1, :] .= 0.0

    #un[:, :, 1] .= 0.0
    #un[:, :, nz-1] .= 0.0
    #vn[:, :, 1] .= 0.0
    #vn[:, :, nz-1] .= 0.0
    #wn[:, :, 1] .= 0.0
    #wn[:, :, nz-1] .= 0.0

    periodic_bc(un, vn, wn, pn, nx, ny, nz, dx, dy, dz, dt, rho)

    #return un, vn, wn

end

#=
function calc_u(u::Array{Float64}, i::Int64, j::Int64, k::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64, F::Int64)

    un = u[i,j,k]-(p[i+1,j,k]-p[i-1,j,k])/(2*dx*rho) + nu*(((u[i+1,j,k]-2*u[i,j,k]+u[i-1,j,k])/(dx^2))+((u[i,j+1,k]-2*u[i,j,k]+u[i,j-1,k])/(dy^2))+((u[i,j,k+1]-2*u[i,j,k]+u[i,j,k-1])/(dz^2))) - u[i,j,k]*((u[i,j,k]-u[i-1,j,k])/(dx)) - v[i,j,k]*((u[i,j,k]-u[i,j-1,k])/(dy)) - w[i,j,k]*((u[i,j,k]-u[i,j,k-1])/(dz))*dt + F*dt

    return un
end

function calc_v(v::Array{Float64}, i::Int64, j::Int64, k::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64)

    vn = v[i,j,k]-(p[i,j+1,k]-p[i,j-1,k])/(2*dy*rho) + nu*(((v[i+1,j,k]-2*v[i,j,k]+v[i-1,j,k])/(dx^2))+(v[i,j+1,k]-2*v[i,j,k]+v[i,j-1,k])/(dy^2))+((v[i,j,k+1]-2*v[i,j,k]+v[i,j,k-1])/(dz^2)) - u[i,j,k]*((v[i,j,k]-v[i-1,j,k])/(dx)) - v[i,j,k]*((v[i,j,k]-v[i,j-1,k])/(dy)) - w[i,j,k]*((v[i,j,k]-v[i,j,k-1])/(dz))*dt

    return vn
end

function calc_w(w::Array{Float64}, i::Int64, j::Int64, k::Int64, rho::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64)

    wn = w[i,j,k]-(p[i,j+1,k]-p[i,j-1,k])/(2*dz*rho) + nu*(((w[i+1,j,k]-2*w[i,j,k]+w[i-1,j,k])/(dx^2))+(w[i,j+1,k]-2*w[i,j,k]+w[i,j-1,k])/(dy^2))+((w[i,j,k+1]-2*w[i,j,k]+w[i,j,k-1])/(dz^2)) - u[i,j,k]*((w[i,j,k]-w[i-1,j,k])/(dx)) - v[i,j,k]*((w[i,j,k]-w[i,j-1,k])/(dy)) - w[i,j,k]*((w[i,j,k]-w[i,j,k-1])/(dz))*dt

    return wn
end

function calc_p(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, dx::Float64, dy::Float64, dz::Float64, rho::Int64, i::Int64, j::Int64, k::Int64)

    p = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) * (((u[i+1,j,k]-u[i-1,j,k])/(2*dx))^2 + ((v[i,j-1,k]-v[i,j+1,k])/(2*dy))^2 + ((w[i,j,k-1]-w[i,j,k-1])/(2*dz))^2 + 2*((u[i,j+1,k]-u[i,j-1,k])/(2*dy))*((v[i+1,j,k]-v[i-1,j,k])/(2*dx)) + 2*((w[i,j+1,k]-w[i,j-1,k])/(2*dy))*((v[i,j,k+1]-v[i,j,k-1])/(2*dz)) + 2*((u[i,j,k+1]-u[i,j,k-1])/(2*dz))*((w[i+1,j,k]-w[i-1,j,k])/(2*dx)))

    return p
end

=#

function periodic_bc(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, p::Array{Float64}, nx::Int64, ny::Int64, nz::Int64, dx::Float64, dy::Float64, dz::Float64, dt::Float64, rho::Int64)

    #Periodic BC at z=0 for u
    u[2:nx-1,2:ny-1,1] = u[2:nx-1,2:ny-1,1]-(p[3:nx,2:ny-1,1]-p[1:nx-2,2:ny-1,1])/(2*dx*rho) + nu*(((u[3:nx,2:ny-1,1]-2*u[2:nx-1,2:ny-1,1]+u[1:nx-2,2:ny-1,1])/(dx^2))+((u[2:nx-1,3:ny,1]-2*u[2:nx-1,3:ny,1]+u[2:nx-1,1:ny-2,1])/(dy^2))+((u[2:nx-1,2:ny-1,2]-2*u[2:nx-1,2:ny-1,1]+u[2:nx-1,2:ny-1, nz-1])/(dz^2))) - u[2:nx-1,2:ny-1,1]*((u[2:nx-1,2:ny-1,1]-u[1:nx-2,2:ny-1,1])/(dx)) - v[2:nx-1,2:ny-1,1]*((u[2:nx-1,2:ny-1,1]-u[2:nx-1,1:ny-2,1])/(dy)) - w[2:nx-1,2:ny-1,1]*((u[2:nx-1,2:ny-1,1]-u[2:nx-1,2:ny-1,nz-1])/(dz))*dt

    #Periodic BC at z=2 for u
    u[2:nx-1,2:ny-1,nz-1] = u[2:nx-1,2:ny-1,nz-1]-(p[3:nx,2:ny-1,nz-1]-p[1:nx-2,2:ny-1,nz-1])/(2*dx*rho) + nu*(((u[3:nx,2:ny-1,nz-1]-2*u[2:nx-1,2:ny-1,nz-1]+u[1:nx-2,2:ny-1,nz-1])/(dx^2))+((u[2:nx-1,3:ny,nz-1]-2*u[2:nx-1,3:ny,nz-1]+u[2:nx-1,1:ny-2,nz-1])/(dy^2))+((u[2:nx-1,2:ny-1, 1]-2*u[2:nx-1,2:ny-1,nz-1]+u[2:nx-1,2:ny-1,nz-2])/(dz^2))) - u[2:nx-1,2:ny-1,nz-1]*((u[2:nx-1,2:ny-1,nz-1]-u[1:nx-2,2:ny-1,nz-1])/(dx)) - v[2:nx-1,2:ny-1,nz-1]*((u[2:nx-1,2:ny-1,nz-1]-u[2:nx-1,1:ny-2,nz-1])/(dy)) - w[2:nx-1,2:ny-1,nz-1]*((u[2:nx-1,2:ny-1,nz-1]-u[2:nx-1,2:ny-1,nz-2])/(dz))*dt
    
    #Periodic BC at z=0 for v
    v[2:nx-1,2:ny-1,1] = v[2:nx-1,2:ny-1,1]-(p[2:nx-1,3:ny,1]-p[2:nx-1,1:ny-2,1])/(2*dy*rho) + nu*(((v[3:nx,2:ny-1,1]-2*v[2:nx-1,2:ny-1,1]+v[1:nx-2,2:ny-1,1])/(dx^2))+(v[2:nx-1,3:ny,1]-2*v[2:nx-1,2:ny-1,1]+v[2:nx-1,1:ny-2,1])/(dy^2))+((v[2:nx-1,2:ny-1,2]-2*v[2:nx-1,2:ny-1,1]+v[2:nx-1,2:ny-1,nz-1])/(dz^2)) - u[2:nx-1,2:ny-1,1]*((v[2:nx-1,2:ny-1,1]-v[1:nx-2,2:ny-1,1])/(dx)) - v[2:nx-1,2:ny-1,1]*((v[2:nx-1,2:ny-1,1]-v[2:nx-1,1:ny-2,1])/(dy)) - w[2:nx-1,2:ny-1,1]*((v[2:nx-1,2:ny-1,1]-v[2:nx-1,2:ny-1,nx-1])/(dz))*dt

    #Periodic BC at z=2 for v
    v[2:nx-1,2:ny-1,nz-1] = v[2:nx-1,2:ny-1,nz-1]-(p[2:nx-1,3:ny,nz-1]-p[2:nx-1,1:ny-2,nz-1])/(2*dy*rho) + nu*(((v[3:nx,2:ny-1,nz-1]-2*v[2:nx-1,2:ny-1,nz-1]+v[1:nx-2,2:ny-1,nz-1])/(dx^2))+(v[2:nx-1,3:ny,nz-1]-2*v[2:nx-1,2:ny-1,nz-1]+v[2:nx-1,1:ny-2,nz-1])/(dy^2))+((v[2:nx-1,2:ny-1,1]-2*v[2:nx-1,2:ny-1,nz-1]+v[2:nx-1,2:ny-1,nz-2])/(dz^2)) - u[2:nx-1,2:ny-1,nz-1]*((v[2:nx-1,2:ny-1,nz-1]-v[1:nx-2,2:ny-1,nz-1])/(dx)) - v[2:nx-1,2:ny-1,nz-1]*((v[2:nx-1,2:ny-1,nz-1]-v[2:nx-1,1:ny-2,nz-1])/(dy)) - w[2:nx-1,2:ny-1,nz-1]*((v[2:nx-1,2:ny-1,nz-1]-v[2:nx-1,2:ny-1,nz-2])/(dz))*dt
    
    #Periodic BC at z=0 for w
    w[2:nx-1,2:ny-1,1] = w[2:nx-1,2:ny-1,1]-(p[2:nx-1,3:ny,1]-p[2:nx-1,1:ny-2,1])/(2*dz*rho) + nu*(((w[3:nx,2:ny-1,1]-2*w[2:nx-1,2:ny-1,1]+w[1:nx-2,2:ny-1,1])/(dx^2))+(w[2:nx-1,3:ny,1]-2*w[2:nx-1,2:ny-1,1]+w[2:nx-1,1:ny-2,1])/(dy^2))+((w[2:nx-1,2:ny-1,2]-2*w[2:nx-1,2:ny-1,1]+w[2:nx-1,2:ny-1,nz-1])/(dz^2)) - u[2:nx-1,2:ny-1,1]*((w[2:nx-1,2:ny-1,1]-w[1:nx-2,2:ny-1,1])/(dx)) - v[2:nx-1,2:ny-1,1]*((w[2:nx-1,2:ny-1,1]-w[2:nx-1,1:ny-2,1])/(dy)) - w[2:nx-1,2:ny-1,1]*((w[2:nx-1,2:ny-1,1]-w[2:nx-1,2:ny-1,nz-1])/(dz))*dt

    #Periodic BC at z=2 for w
    w[2:nx-1,2:ny-1,nz-1] = w[2:nx-1,2:ny-1,nz-1]-(p[2:nx-1,3:ny,nz-1]-p[2:nx-1,1:ny-2,nz-1])/(2*dz*rho) + nu*(((w[3:nx,2:ny-1,nz-1]-2*w[2:nx-1,2:ny-1,nz-1]+w[1:nx-2,2:ny-1,nz-1])/(dx^2))+(w[2:nx-1,3:ny,nz-1]-2*w[2:nx-1,2:ny-1,nz-1]+w[2:nx-1,1:ny-2,nz-1])/(dy^2))+((w[2:nx-1,2:ny-1,1]-2*w[2:nx-1,2:ny-1,nz-1]+w[2:nx-1,2:ny-1,nz-2])/(dz^2)) - u[2:nx-1,2:ny-1,nz-1]*((w[2:nx-1,2:ny-1,nz-1]-w[1:nx-2,2:ny-1,nz-1])/(dx)) - v[2:nx-1,2:ny-1,nz-1]*((w[2:nx-1,2:ny-1,nz-1]-w[2:nx-1,1:ny-2,nz-1])/(dy)) - w[2:nx-1,2:ny-1,nz-1]*((w[2:nx-1,2:ny-1,nz-1]-w[2:nx-1,2:ny-1,nz-2])/(dz))*dt

    #Periodic BC at z=0 for p
    p[2:nx-1,2:ny-1,1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) * (((u[3:nx,2:ny-1,1]-u[1:nx-2,2:ny-1,1])/(2*dx))^2 + ((v[2:nx-1,1:ny-2,1]-v[2:nx-1,3:ny,1])/(2*dy))^2 + ((w[2:nx-1,2:ny-1,nz-1]-w[2:nx-1,2:ny-1,nz-1])/(2*dz))^2 + 2*((u[2:nx-1,3:ny,1]-u[2:nx-1,1:ny-2,1])/(2*dy))*((v[3:nx,2:ny-1,1]-v[1:nx-2,2:ny-1,1])/(2*dx)) + 2*((w[2:nx-1,3:ny,1]-w[2:nx-1,1:ny-2,1])/(2*dy))*((v[2:nx-1,2:ny-1,2]-v[2:nx-1,2:ny-1,nz-1])/(2*dz)) + 2*((u[2:nx-1,2:ny-1,2]-u[2:nx-1,2:ny-1,nz-1])/(2*dz))*((w[3:nx,2:ny-1,1]-w[1:nx-2,2:ny-1,1])/(2*dx)))
    
    #Periodic BC at z=2 for p
    p[2:nx-1,2:ny-1,nz-1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) * (((u[3:nx,2:ny-1,nz-1]-u[1:nx-2,2:ny-1,nz-1])/(2*dx))^2 + ((v[2:nx-1,1:ny-2,nz-1]-v[2:nx-1,3:ny,nz-1])/(2*dy))^2 + ((w[2:nx-1,2:ny-1,nz-2]-w[2:nx-1,2:ny-1,nz-2])/(2*dz))^2 + 2*((u[2:nx-1,3:ny,nz-1]-u[2:nx-1,1:ny-2,nz-1])/(2*dy))*((v[3:nx,2:ny-1,nz-1]-v[1:nx-2,2:ny-1,nz-1])/(2*dx)) + 2*((w[2:nx-1,3:ny,nz-1]-w[2:nx-1,1:ny-2,nz-1])/(2*dy))*((v[2:nx-1,2:ny-1,1]-v[2:nx-1,2:ny-1,nz-2])/(2*dz)) + 2*((u[2:nx-1,2:ny-1,1]-u[2:nx-1,2:ny-1,nz-2])/(2*dz))*((w[3:nx,2:ny-1,nz-1]-w[1:nx-2,2:ny-1,nz-1])/(2*dx)))



end

function meshgrid3d(xin,yin,zin)
    nx=length(xin)
    ny=length(yin)
    nz=length(zin)
    xout=zeros(nz,ny,nx)
    yout=zeros(nz,ny,nx)
    zout=zeros(nz,ny,nz)
    
    
    for kx=1:nx
        for jx=1:ny
            for ix=1:nz
                xout[ix,jx,kx]=xin[kx]
                yout[ix,jx,kx]=yin[jx]
                zout[ix,jx,kx]=zin[ix]
            end
        end
    end
    return (x=xout, y=yout, z=zout)
end

function map3d(u, v, w, X, Y, Z)
    
    uvw = [Point3f(u[1,1,1], v[1,1,1], w[1,1,1])]
    xyz = [Point3f(X[1,1,1], Y[1,1,1], Z[1,1,1])]
    strength = vec(sqrt.(u[2:nx-1, 2:ny-1, 2:nz-1] .^ 2 .+ v[2:nx-1, 2:nx-1, 2:nz-1] .^ 2 .+ w[2:nx-1, 2:ny-1, 2:nz-1]))

    for i in range(1, stop=nx)
        for j in range(1, stop=ny)
            for k in range(1, stop=nz)
                if (i != 1 && j != 1 && k != 1)
                    push!(uvw, Point3f(u[i,j,k], v[i,j,k], w[i,j,k]))
                    push!(xyz, Point3f(X[i,j,k], Y[i,j,k], Z[i,j,k]))
                end
            end
        end
    end

    return uvw, xyz, strength

end

const nx = 41
const ny = 41
const nz = 41
const nt = 10
nit = 2
const c = 1
const dx = 2 / (nx - 1)
const dy = 2 / (ny - 1)
const dz = 2 / (nz - 1)

x = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = Array{Float64}([i for i in range(0.0, stop=2.0, step=dy)])
z = Array{Float64}([i for i in range(0.0, stop=2.0, step=dz)])

const rho = 1
const nu = 0.1
const dt = 0.000001
const F = 1

u = zeros(nx, ny, nz)
v = zeros(nx, ny, nz)
w = zeros(nx, ny, nz)
p = zeros(nx, ny, nz)
un = zeros(nx, ny, nz)
vn = zeros(nx, ny, nz)
wn = zeros(nx, ny, nz)
pn = zeros(nx, ny, nz)

X,Y,Z = meshgrid3d(x,y,z)

n = 0

#f = Figure(resolution = (800, 800))
#Axis(f[1,1], backgroundcolor = "white")

#record(f, "Channel Flow.gif", 1:200) do l
while n<nit 
    run(u, v, w, p, un, vn, wn, pn, nx, ny, nz, rho, dx, dy, dz, dt, F, X, Y, Z)
    #print(uvw)
    #arrows!(xyz, uvw)
    global n = n + 1
end
uvw, xyz, strength = map3d(un, vn, wn, X, Y, Z)
arrows(xyz, uvw, lengtscale = 0.001, arrowsize=0.03)