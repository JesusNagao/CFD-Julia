
using Plots
using Base.Threads
using DelimitedFiles

mutable struct Wave2D
    nx::Int64
    ny::Int64
    nt::Int64
    u::Array{Float64,3} 
    v::Array{Float64,3} 
end

Wave2D(nx,ny,nt) = Wave2D(nx,ny,nt,ones(nx,ny,nt),ones(nx,ny,nt)) 

function Convection(w::Wave2D, exportcsv = false)
    nx::Int64 = w.nx
    ny::Int64 = w.ny
    nt::Int64 = w.nt
    c::Float64 = 1
    σ::Float64 = 0.2
    dx::Float64 = 2/(nx-1)
    dy::Float64 = 2/(ny-1) 
    dt::Float64 = σ*dx

    X = LinRange(0,2,nx)
    Y = LinRange(0,2,nx)

    #Condiciones iniciales
    w.u[floor(Int64,0.5/dx):floor(Int64,1/ dx + 1),floor(Int64,0.5/dy):floor(Int64,1/ dy + 1),1] .= 2
    w.v[floor(Int64,0.5/dx):floor(Int64,1/ dx + 1),floor(Int64,0.5/dy):floor(Int64,1/ dy + 1),1] .= 2
    
    convect = @animate for n in range(1,stop=nt-1)

        @threads for i in range(2,stop = nx)

            @threads for j in range(2,stop =ny)

                w.u[i,j,n+1] = w.u[i,j,n] - w.u[i,j,n]*dt/dx*(w.u[i,j,n]-w.u[i-1,j,n]) - w.v[i,j,n]*dt/dy*(w.u[i,j,n]-w.u[i,j-1,n])

                w.v[i,j,n+1] = w.v[i,j,n] - w.u[i,j,n]*dt/dx*(w.v[i,j,n]-w.v[i-1,j,n]) - w.v[i,j,n]*dt/dy*(w.v[i,j,n]-w.v[i,j-1,n])

            end
        end
        
        w.u[1, :,n] .= 1
        w.u[end, :,n] .= 1
        w.u[:, 1,n] .= 1
        w.u[:, end,n] .= 1

        plot(X,Y,w.u[:,:,n],st = :surface, xlims = (0,2), ylims = (0,2), zlims = (0,2), camera=(45,70), c = cgrad(:ice))

    end
    
    if exportcsv == true
        writedlm("conv.csv", w.u, ',')
    end

    gif(convect, "conv.gif", fps = 1000)

end

conv1 = Wave2D(101,101,160)
@time Convection(conv1)