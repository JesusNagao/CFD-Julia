include("init_speed.jl")
include("time.jl")
#using Plots;

const nx = 81;
const ny = 81;
const nt = 100;
const c = 1;
const dx = 2 / (nx-1);
const dy = 2 / (ny-1);
const sigma = 0.2;
const dt = sigma * dx;

x = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])
y = Array{Float64}([i for i in range(0.0, stop=2.0, step=dx)])

u = Matrix{Float64}(undef, nx, ny)
un = Matrix{Float64}(undef, nx, ny)

init_speed(u, dx, dy)
time_iteration(u, un, nt, c, dt, dx, dy)
#surface(x, y, u)

