using GLMakie
using CUDA
include("GPU Operations.jl")

const nx = 41
const ny = 41
const nt = 10
const nit = 51
const c = 1
const dx = 2 / (nx - 1)
const dy = 2 / (ny - 1)

x = CuArray{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = CuArray{Float64}([i for i in range(0.0, stop=2.0, step=dy)])

const rho = 1
const nu = 0.1
const F = 1
const dt = 0.01

u = CuArray{Float64}(undef, nx, ny)
u = CUDA.zeros(nx, ny)
v = CuArray{Float64}(undef, nx, ny)
v = CUDA.zeros(nx, ny)
p = CuArray{Float64}(undef, nx, ny)
p = CUDA.zeros(nx, ny)
b = CuArray{Float64}(undef, nx, ny)
b = CUDA.zeros(nx, ny)
un = CuArray{Float64}(undef, nx, ny)
un = CUDA.zeros(nx, ny)
vn = CuArray{Float64}(undef, nx, ny)
vn = CUDA.zeros(nx, ny)
pn = CuArray{Float64}(undef, nx, ny)
pn = CUDA.zeros(nx, ny)

run_gpu()