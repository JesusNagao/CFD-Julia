const nx = 31
const ny = 31
const nt = 17
const nu = 0.05
const dx = 2/(nx-1)
const dy = 2/(ny-1)
const sigma = 0.25
const dt = sigma * dx * dy / nu

x = Array{Float64}([i for i in range(0, stop=2, step=dx)])
y = Array{Float64}([i for i in range(0, stop=2, step=dy)])