using GLMakie,LinearAlgebra
include("3DNavier.jl")

#Produces a 3D vector plot of a given field

#Solving NS equations on given conditions
nx = 41
ny = 41
nz = 41
nt = 500
nit = 50

c = 1

dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
dz = 2 / (nz - 1)

x = LinRange(0, 2, nx)
y = LinRange(0, 2, ny)
z = LinRange(0, 2, nz)

rho = 1.
nu = .1
dt = .001

u = zeros((nx, ny, nz))
v = zeros((nx, ny, nz))
w = zeros((nx, ny, nz))

p = zeros((nx, ny, nz)) 
b = zeros((nx, ny, nz))

X,Y,Z = meshgrid(x,y,z)

u, v, w, p = cavity_flow(nt, u, v, w, dt, dx, dy, dz, p, rho, nu)

#Correr una vez obtenidos los vectores u,v,p

#fig = CairoMakie.Figure(resolution = (800, 600))

function map3d(u, v, w, X, Y, Z)
    
    uvw = [Point3f0(u[1,1,1], v[1,1,1], w[1,1,1])]
    xyz = [Point3f0(X[1,1,1], Y[1,1,1], Z[1,1,1])]

    for i in range(1, stop=nx)
        for j in range(1, stop=ny)
            for k in range(1, stop=nz)
                if (i != 1 && j != 1 && k != 1)
                    push!(uvw, Point3f0(u[i,j,k], v[i,j,k], w[i,j,k]))
                    push!(xyz, Point3f0(X[i,j,k], Y[i,j,k], Z[i,j,k]))
                end
            end
        end
    end

    return uvw, xyz

end

#CairoMakie.contourf(x, y, p[:,:,20], levels = 5, linewidth = 3)
#CairoMakie.arrows(x,y,u[:,:,50],v[:,:,50], arrowcolor = leng)

uvw, xyz= map3d(u, v, w, X, Y, Z)

#Normalizing
uvwnorm = [norm(vec) != 0 ? vec/norm(vec) : vec for vec in uvw]

arrows(xyz[1:3:end], uvw[1:3:end], lengthscale = 0.1, arrowsize=0.01,linewidth = 0.05,linecolor=norm.(uvw[1:3:end]))

#CairoMakie.save("N2D.png", fig)