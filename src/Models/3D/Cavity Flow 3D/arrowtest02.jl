using GLMakie, Zygote, LinearAlgebra
x = y = -2:0.1:2
g(x, y) = 4*sin(x)*exp(-(x^2+y^2));
z = g.(x, y')

fig = Figure(resolution=(1500,1500))
ax = Axis3(fig[1,1])
surface!(x, y, g.(x, y'))
pos = vec(Point3.(x, y', g.(x, y')))
grad = vec((p->-Vec3(p..., 0.) / 8).(gradient.(g, x, y')))
#cmap = RGBAf0.(to_colormap(:viridis), 0.6)

arrows!(vec(Point3.(x, y', g.(x, y'))), grad;
        arrowsize = 0.05,
        linecolor=norm.(grad),
        linewidth = 0.01)#,
        #colormap=cmap)
fig
