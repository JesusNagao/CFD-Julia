"""
General square grid in N dimensions along with a velocity field
"""
struct SquareGrid{N}
    𝐱::Array{Float64,N} 
    𝐮::Array{Float64,N} 
end