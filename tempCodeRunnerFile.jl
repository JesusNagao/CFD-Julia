nx = 41
ny = 41
nz = 41
nt = 1500
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