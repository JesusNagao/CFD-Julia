using CairoMakie
using LinearAlgebra

#Correr una vez obtenidos los vectores u,v,p

#fig = CairoMakie.Figure(resolution = (800, 600))

x = LinRange(0,2,41)
y = LinRange(0,2,41)
leng = vec(norm.(Vec2f0.(u[:,:,1],v[:,:,1])))

#CairoMakie.contourf(x, y, p[:,:,20], levels = 5, linewidth = 3)
CairoMakie.arrows(x,y,u[:,:,50],v[:,:,50], arrowcolor = leng)

#CairoMakie.save("N2D.png", fig)
