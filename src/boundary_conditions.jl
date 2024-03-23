"""Periodic and static boundary condition definitions for 2 and 3 dimensional grids"""

"""Apply periodic boundary to a grid."""
function periodic_boundary(grid::SquareGrid{2})
    """Apply periodic boundary to a 2D grid."""
    grid.𝐮[1, :] = grid.𝐮[end-1, :]
    grid.𝐮[:, 1] = grid.𝐮[:, end-1]
    grid.𝐯[1, :] = grid.𝐯[end-1, :]
    grid.𝐯[:, 1] = grid.𝐯[:, end-1]
    
end

function periodic_boundary(grid::SquareGrid{3})
    """Apply periodic boundary to a grid."""
    grid.𝐮[1, :, :] = grid.𝐮[end-1, :,:]
    grid.𝐮[:, 1, :] = grid.𝐮[:, end-1, :]
    grid.𝐯[1, :, :] = grid.𝐯[end-1, :, :]
    grid.𝐯[end, :, :] = grid.𝐯[2, :, :]
    grid.𝐰[1, :, :] = grid.𝐰[end-1, :, :]
    grid.𝐰[end, :, :] = grid.𝐰[2, :, :]

end

"""Apply static boundary conditions to a grid"""

function static_boundary(grid::SquareGrid{2})
    """Apply static boundary conditions to a 2D grid"""
    grid.𝐮[1, :] .= 0
    grid.𝐮[end, :] .= 0
    grid.𝐮[:, 1] .= 0
    grid.𝐮[:, end] .= 0
    grid.𝐯[1, :] .= 0
    grid.𝐯[end, :] .= 0
    grid.𝐯[:, 1] .= 0
    grid.𝐯[:, end] .= 0

end

function static_boundary(grid::SquareGrid{3})
    """Apply static boundary conditions to a 3D grid"""
    grid.𝐮[1, :] .= 0
    grid.𝐮[end, :] .= 0
    grid.𝐮[:, 1] .= 0
    grid.𝐮[:, end] .= 0
    grid.𝐯[1, :] .= 0
    grid.𝐯[end, :] .= 0
    grid.𝐯[:, 1] .= 0
    grid.𝐯[:, end] .= 0
    grid.𝐰[1, :] .= 0
    grid.𝐰[end, :] .= 0
    grid.𝐰[:, 1] .= 0
    grid.𝐰[:, end] .= 0
    
end