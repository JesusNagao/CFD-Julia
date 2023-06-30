using LinearAlgebra
using BenchmarkTools
using GLMakie
using Makie
include("3DNavier.jl")
include("3DNavierthreads.jl")

#Produces a 3D vector plot of a given field

#Solving NS equations on given conditions
nx = 41
ny = 41
nz = 41
nt = 700
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

#Plotting with Makie 

#Correr una vez obtenidos los vectores u,v,w, p

u,v,w,p = cavity_flow(nt, u, v, w, dt, dx, dy, dz, p, rho, nu)

#leng = vec(norm.(Vec3f0.(u,v,w)))

arrows(x,y,z,u,v,w)

#CairoMakie.save("N2D.png", fig)

#@btime cavity_flow(nt, u, v, w, dt, dx, dy, dz, p, rho, nu)
#@btime cavity_flow_threads(nt, u, v, w, dt, dx, dy, dz, p, rho, nu)