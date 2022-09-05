using PyPlot
using LinearAlgebra
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

@time u, v, w, p = cavity_flow(nt, u, v, w, dt, dx, dy, dz, p, rho, nu)

#Plotting the obtained field

fig = figure()
ax = fig.gca(projection="3d")

ax.quiver(x,y,z, u,v,w, length = 0.5)
ax.view_init(0, 45)

ax.set_title("Solución a ecuaciones de momentum")
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

fig