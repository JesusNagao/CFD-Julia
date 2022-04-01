using Plots

function init(u::Matrix{Float64}, v::Matrix{Float64}, nx::Int64, ny::Int64)

    for i in range(1, stop=nx)
        for j in range(1, stop=ny)

            for i in range(1, stop=size(u)[1])
                for j in range(1, stop=size(u)[2])
                    if (x[i] < 1 && x[i] > 0.5 && y[j] < 1 && y[j] > 0.5 )
        
                        u[i, j] = 2
                        v[i, j] = 2
        
                    end
                end
            end

        end
    end

end

function animate(u::Matrix{Float64}, v::Matrix{Float64}, un::Matrix{Float64}, vn::Matrix{Float64}, nx::Int64, ny::Int64, dx::Float64, dy::Float64, dt::Float64, nu::Float64)

    @gif for k in range(1, stop=nt)

        un = u
        vn = v

        burger(u, v, un, vn, nx, ny, dx, dy, dt, nu)
        surface(x, y, u, xlim=(0,2), ylim=(0,2), zlim=(0,2), legend=false)

    end

end

function burger(u::Matrix{Float64}, v::Matrix{Float64}, un::Matrix{Float64}, vn::Matrix{Float64}, nx::Int64, ny::Int64, dx::Float64, dy::Float64, dt::Float64, nu::Float64)

    Threads.@threads for i  in range(2, stop=nx-1)
        Threads.@threads for j in range(2, stop=ny-1)
            u[i, j] = un[i,j] - dt*un[i,j]*(un[i,j]-un[i-1, j])/dx - dt*vn[i,j]*(un[i,j]-un[i, j-1])/dy + nu*dt*(un[i+1, j]-2*un[i,j]+un[i-1,j])/(dx^2) + nu*dt*(un[i,j+1]-2*un[i,j]+un[i,j-1])/(dy^2)
            v[i, j] = vn[i,j] - dt*un[i,j]*(vn[i,j]-vn[i-1, j])/dx - dt*vn[i,j]*(vn[i,j]-vn[i, j-1])/dy + nu*dt*(vn[i+1, j]-2*vn[i,j]+vn[i-1,j])/(dx^2) + nu*dt*(vn[i,j+1]-2*vn[i,j]+vn[i,j-1])/(dy^2)
            
            #v[i, j] = vn[i,j] + nu*dt/(dx^2)*(vn[i+1, j] - 2*vn[i,j] + vn[i-1, j]) + nu*dt/(dy^2)*(vn[i, j+1]-2*vn[i,j]+vn[i,j-1])
        end
    end

end