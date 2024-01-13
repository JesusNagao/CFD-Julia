include("Operations.jl")

const nx = 41
const ny = 41
const nt = 120
const c = 1
const dx = 2/(nx-1)
const dy = 2/(ny-1)
const sigma = 0.0009
const nu = 0.01
const dt = sigma * dx * dy / nu

x = Array{Float64}([i for i in range(0, stop=2, step=dx)])
y = Array{Float64}([i for i in range(0, stop=2, step=dy)])

u = Matrix{Float64}(undef, nx, ny)
u = zeros(nx, ny)
v = Matrix{Float64}(undef, nx, ny)
v = zeros(nx, ny)

un = Matrix{Float64}(undef, nx, ny)
un = zeros(nx, ny)
vn = Matrix{Float64}(undef, nx, ny)
vn = zeros(nx, ny)
comb = Matrix{Float64}(undef, nx, ny)
comb = zeros(nx, ny)

init(u, v, nx, ny)
animate(u, v, un, vn, nx, ny, dx, dy, dt, nu)