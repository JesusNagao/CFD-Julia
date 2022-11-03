using GLMakie
include("Operations.jl")
#using LinearAlgebra

const nx = 41
const ny = 41
const nz = 41
const nt = 10
const c = 1
const dx = 2 / (nx - 1)
const dy = 2 / (ny - 1)
const dz = 2 / (nz - 1)
const rho = 1
const nu = 0.1
const F = 1

nit = 10
dt = 0.001

x = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = Array{Float64}([i for i in range(0.0, stop=2.0, step=dy)])
z = Array{Float64}([i for i in range(0.0, stop=2.0, step=dz)])

X,Y,Z = meshgrid3d(x,y,z)

u = zeros(nx, ny, nz)
v = zeros(nx, ny, nz)
w = zeros(nx, ny, nz)
p = zeros(nx, ny, nz)
un = zeros(nx, ny, nz)
vn = zeros(nx, ny, nz)
wn = zeros(nx, ny, nz)
pn = zeros(nx, ny, nz)

f = Figure(resolution = (800, 800))

initiate(x, y, z, un, vn, wn, nx, ny, nz)
record(f, "Box.gif", 1:50) do i
    #print("Hola")
    Axis3(f[1,1], backgroundcolor = "black")
    iterate(u,v,w,un,vn,wn,p,nx,ny, nz,nu,dx,dy,dz,dt,F)
    uvw, xyz, strength = map3d(un, vn, wn, X, Y, Z)
    arrows!(xyz, uvw, lengthscale = 0.001, arrowcolor = strength, linecolor = strength)
end

#iterate(u,v,w,un,vn,wn,p,nx,ny, nz,nu,dx,dy,dz,dt,F)

