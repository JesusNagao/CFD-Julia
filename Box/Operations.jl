using GLMakie

#module Ops

#export iterate, initiate, map3d, meshgrid3d

#using GLMakie, LinearAlgebra

function map3d(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, X::Array{Float64}, Y::Array{Float64}, Z::Array{Float64})
    
    uvw = [Point3f(u[1,1,1], v[1,1,1], w[1,1,1])]
    xyz = [Point3f(X[1,1,1], Y[1,1,1], Z[1,1,1])]
    #strength = norm.(vec(sqrt.(u[2:nx-1, 2:ny-1, 2:nz-1] .^ 2 .+ v[2:nx-1, 2:nx-1, 2:nz-1] .^ 2 .+ w[2:nx-1, 2:ny-1, 2:nz-1])))
    strength = [norm.(sqrt.(u[1, 1, 1] .^ 2 .+ v[1, 1, 1] .^ 2 .+ w[1, 1, 1]))]

    for i in range(1, stop=nx)
        for j in range(1, stop=ny)
            for k in range(1, stop=nz)
                if (i != 1 && j != 1 && k != 1)
                    push!(uvw, Point3f(u[i,j,k], v[i,j,k], w[i,j,k]))
                    push!(xyz, Point3f(X[i,j,k], Y[i,j,k], Z[i,j,k]))
                    push!(strength, norm.(sqrt.(u[i, j, k] .^ 2 .+ v[i, j, k] .^ 2 .+ w[i, j, k])))
                end
            end
        end
    end

    return uvw, xyz, strength

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


function iterate(u::Array{Float64}, v::Array{Float64}, w::Array{Float64}, un::Array{Float64}, vn::Array{Float64}, wn::Array{Float64}, p::Array{Float64}, nx::Int64, ny::Int64, nz::Int64, nu::Float64, dx::Float64, dy::Float64, dz::Float64, dt::Float64, F::Int64)
    
    pressure(p, un, vn, wn, nx, ny, nz, dx, dy, dz)
    
    u[:,:,:] = un[:,:,:]
    v[:,:,:] = vn[:,:,:]
    w[:,:,:] = wn[:,:,:]
    #p[:,:,:] = pn[:,:,:]

    #pn[2:nx-1, 2:ny-1, 2:nz-1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) * (((u[3:nx,2:ny-1,2:nz-1]-u[1:nx-2,2:ny-1,2:nz-1])/(2*dx)).^2 .+ ((v[2:nx-1,1:ny-2,2:nz-1]-v[2:nx-1,3:ny,2:nz-1])/(2*dy)).^2 .+ ((w[2:nx-1,2:ny-1,1:nz-2]-w[2:nx-1,2:ny-1,1:nz-2])/(2*dz)).^2 .+ 2*((u[2:nx-1,3:ny,2:nz-1]-u[2:nx-1,1:ny-2,2:nz-1])/(2*dy)).*((v[3:nx,2:ny-1,2:nz-1]-v[1:nx-2,2:ny-1,2:nz-1])/(2*dx)) + 2*((w[2:nx-1,3:ny,2:nz-1]-w[2:nx-1,1:ny-2,2:nz-1])/(2*dy)).*((v[2:nx-1,2:ny-1,3:nz]-v[2:nx-1,2:ny-1,1:nz-2])/(2*dz)) + 2*((u[2:nx-1,2:ny-1,3:nz]-u[2:nx-1,2:ny-1,1:nz-2])/(2*dz)).*((w[3:nx,2:ny-1,2:nz-1]-w[1:nx-2,2:ny-1,2:nz-1])/(2*dx)))
    un[2:nx-1, 2:ny-1, 2:nz-1] = u[2:nx-1, 2:ny-1, 2:nz-1] .+ (-(p[3:nx, 2:ny-1, 2:nz-1] - p[1:nx-2, 2:ny-1, 2:nz-1])/(2*rho*dx) + nu*((u[3:nx, 2:ny-1, 2:nz-1] - 2*u[2:nx-1, 2:ny-1, 2:nz-1] + u[1:nx-2, 2:ny-1, 2:nz-1])/(dx^2) + (u[2:nx-1, 3:ny, 2:nz-1] - 2*u[2:nx-1, 2:ny-1, 2:nz-1] + u[2:nx-1, 1:ny-2, 2:nz-1])/(dy^2) + (u[2:nx-1, 2:ny-1, 3:nx] - 2*u[2:nx-1, 2:ny-1, 2:nz-1] + u[2:nx-1, 2:ny-1, 1:nz-2])/(dz^2)) - u[2:nx-1, 2:ny-1, 2:nz-1].*(u[2:nx-1, 2:ny-1, 2:nz-1] - u[1:nx-2, 2:ny-1, 2:nz-1])/(dx) - v[2:nx-1, 2:ny-1, 2:nz-1].*(u[2:nx-1, 2:ny-1, 2:nz-1] - u[2:nx-1, 1:ny-2, 2:nz-1])/(dy) - w[2:nx-1, 2:ny-1, 2:nz-1].*(u[2:nx-1, 2:ny-1, 2:nz-1] - u[2:nx-1, 2:ny-1, 1:nz-2])/(dz))*dt .+ F*dt
    vn[2:nx-1, 2:ny-1, 2:nz-1] = v[2:nx-1, 2:ny-1, 2:nz-1] .+ (-(p[2:nx-1, 3:ny, 2:nz-1] - p[2:nx-1, 1:ny-2, 2:nx-1])/(2*dy*rho) + nu*((v[3:nx, 2:ny-1, 2:nz-1] - 2*v[2:nx-1, 2:ny-1, 2:nz-1] + v[1:nx-2, 2:ny-1, 2:nz-1])/(dx^2) + (v[2:nx-1, 3:ny, 2:nz-1] - 2*v[2:nx-1, 2:ny-1, 2:nz-1] + v[2:nx-1, 1:ny-2, 2:nz-1])/(dy^2) + (v[2:nx-1, 2:ny-1, 3:nz] - 2*v[2:nx-1, 2:ny-1, 2:nz-1] + v[2:nx-1, 2:ny-1, 1:nz-2])/(dz^2)) - u[2:nx-1, 2:ny-1, 2:nz-1].*(v[2:nx-1, 2:ny-1, 2:nz-1] - v[1:nx-2, 2:ny-1, 2:nz-1])/(dx) - v[2:nx-1, 2:ny-1, 2:nz-1].*(v[2:nx-1, 2:ny-1, 2:nz-1] - v[2:nx-1, 1:ny-2, 2:nz-1])/(dy) - w[2:nx-1, 2:ny-1, 2:nz-1].*(v[2:nx-1, 2:ny-1, 2:nz-1] - v[2:nx-1, 2:ny-1, 1:nz-2])/(dz))*dt
    wn[2:nx-1, 2:ny-1, 2:nz-1] = w[2:nx-1, 2:ny-1, 2:nz-1] .+ (-(p[2:nx-1, 2:ny-1, 3:nz] - p[2:nx-1, 2:ny-1, 1:nz-2])/(2*dz*rho) + nu*((w[3:nx, 2:ny-1, 2:nz-1] - 2*w[2:nx-1, 2:ny-1, 2:nz-1] + w[1:nx-2, 2:ny-1, 2:nz-1])/(dx^2) + (w[2:nx-1, 3:ny, 2:nz-1] - 2*w[2:nx-1, 2:ny-1, 2:nz-1] + w[2:nx-1, 1:ny-2, 2:nz-1])/(dy^2) + (w[2:nx-1, 2:ny-1, 3:nz] - 2*w[2:nx-1, 2:ny-1, 2:nz-1] + w[2:nx-1, 2:ny-1, 1:nz-2])/(dz^2)) - u[2:nx-1, 2:ny-1, 2:nz-1].*(w[2:nx-1, 2:ny-1, 2:nz-1] - w[1:nx-2, 2:ny-1, 2:nz-1])/(dx) - v[2:nx-1, 2:ny-1, 2:nz-1].*(w[2:nx-1, 2:ny-1, 2:nz-1] - w[2:nx-1, 1:ny-2, 2:nz-1])/(dy) - w[2:nx-1, 2:ny-1, 2:nz-1].*(w[2:nx-1, 2:ny-1, 2:nz-1] - w[2:nx-1, 2:ny-1, 1:nz-2])/(dz))*dt

    boundaries(un, vn, wn, nx, ny, nz)

    #print(un)

end

function pressure(p, u, v, w, nx, ny, nz, dx, dy, dz)
    p[nx-1, :,:] = p[nx-2, :,:]
    p[1, :, :] = p[2, :, :]
    p[:, ny-1,:] = p[:, ny-2,:]
    p[:, 1, :] = p[:, 2, :]
    p[:, :, nz-1] = p[:, :, nz-2]
    p[:, :, 1] = p[:, :, 2]

    p[2:nx-1, 2:ny-1, 2:nz-1] = ((dx^2*dy^2*dz^2*rho)/(2*(dy^2*dz^2+dx^2*dz^2+dy^2*dx^2))) * (((u[3:nx,2:ny-1,2:nz-1]-u[1:nx-2,2:ny-1,2:nz-1])/(2*dx)).^2 .+ ((v[2:nx-1,1:ny-2,2:nz-1]-v[2:nx-1,3:ny,2:nz-1])/(2*dy)).^2 .+ ((w[2:nx-1,2:ny-1,1:nz-2]-w[2:nx-1,2:ny-1,1:nz-2])/(2*dz)).^2 .+ 2*((u[2:nx-1,3:ny,2:nz-1]-u[2:nx-1,1:ny-2,2:nz-1])/(2*dy)).*((v[3:nx,2:ny-1,2:nz-1]-v[1:nx-2,2:ny-1,2:nz-1])/(2*dx)) + 2*((w[2:nx-1,3:ny,2:nz-1]-w[2:nx-1,1:ny-2,2:nz-1])/(2*dy)).*((v[2:nx-1,2:ny-1,3:nz]-v[2:nx-1,2:ny-1,1:nz-2])/(2*dz)) + 2*((u[2:nx-1,2:ny-1,3:nz]-u[2:nx-1,2:ny-1,1:nz-2])/(2*dz)).*((w[3:nx,2:ny-1,2:nz-1]-w[1:nx-2,2:ny-1,2:nz-1])/(2*dx)))

end


function boundaries(un, vn, wn, nx, ny, nz)
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

    un[:, :, 1] .= 0.0
    un[:, :, nz-1] .= 0.0
    vn[:, :, 1] .= 0.0
    vn[:, :, nz-1] .= 0.0
    wn[:, :, 1] .= 0.0
    wn[:, :, nz-1] .= 0.0
end

function initiate(x,y,z,u,v,w,nx,ny,nz)

    for i in range(1, stop=nx-1)
        for j in range(1, stop=ny-1)
            for k in range(1, stop=nz-1)
                u[i,j,k] = exp((-((x[i]-1)^2)-((y[j]-1)^2)-((z[k]-1)^2))/sqrt(2))
                v[i,j,k] = exp((-((x[i]-1)^2)-((y[j]-1)^2)-((z[k]-1)^2))/sqrt(2))
                w[i,j,k] = exp((-((x[i]-1)^2)-((y[j]-1)^2)-((z[k]-1)^2))/sqrt(2))
            end
        end
    end
    
end

#end