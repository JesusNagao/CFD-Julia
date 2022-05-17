using CairoMakie
using LinearAlgebra
using Colors

#Correr una vez obtenidos los vectores u,v,p

#fig = CairoMakie.Figure(resolution = (800, 600))

x = LinRange(0,2,41)
y = LinRange(0,2,41)
leng = vec(norm.(Vec2f0.(u[:,:,1],v[:,:,1])))

fig, ax, cf = contour(x, y, p[:,:,1])
ar = arrows!(x,y,u[:,:,1],v[:,:,1],lengthscale = 0.1, colormap = :blues)

record(fig, "CavityFlow.mp4", 1:length(p[1,1,:])) do i

    #leng = vec(norm.(Vec2f0.(u[:,:,i],v[:,:,i])))

    cf[1] = x 
    cf[2] = y 
    cf[3] =  p[:,:,i] 

    ar[1] = x
    ar[2] = y
    ar[3] = u[:,:,i]
    ar[4] = v[:,:,i]

    #xlims!(ax,[0,2]) 
    #ylims!(ax,[0,2])

end

#CairoMakie.save("N2D.png", fig)
