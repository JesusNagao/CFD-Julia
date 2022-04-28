using GLMakie
include("Operations.jl")

const nx = 41
const ny = 41
const nt = 10
const nit = 51
const c = 1
const dx = 2 / (nx - 1)
const dy = 2 / (ny - 1)

x = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = Array{Float64}([i for i in range(0.0, stop=2.0, step=dy)])

const rho = 1
const nu = 0.1
const F = 1
const dt = 0.01

u = Matrix{Float64}(undef, nx, ny)
u = zeros(nx, ny)
v = Matrix{Float64}(undef, nx, ny)
v = zeros(nx, ny)
p = Matrix{Float64}(undef, nx, ny)
p = zeros(nx, ny)
b = Matrix{Float64}(undef, nx, ny)
b = zeros(nx, ny)
un = Matrix{Float64}(undef, nx, ny)
un = zeros(nx, ny)
vn = Matrix{Float64}(undef, nx, ny)
vn = zeros(nx, ny)


