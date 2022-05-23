module Plastic
using Reexport, Printf

# 调用者一定也会用到这些包
@reexport using Ferrite
using Tensors
@reexport using WriteVTK

getdim = Ferrite.getdim

# 调用者不一定会用到这些包
using LinearAlgebra
using SparseArrays

import PointInPoly: pinpoly

export  
    J2Plasticity,
    vonMises,
    PlasticStructure,
    BoundaryPoly,
    NewtonRaphsonSolver,
    NewmarkSolver,
    ExplicitSolver
    
export
    generate_grid,
    add_bc!,
    advance!,
    fetch_poly,
    fetch_data,
    save_to_vtk,
    time_step!
    

# 不知道有什么用
# using Base: @propagate_inbounds

# basic types
include("base.jl")

# helper functions
include("geometry.jl")
include("assembler_static.jl")
include("assembler_dynamic.jl")
include("creater.jl")
include("solver.jl")
include("fetchers.jl")
include("postprocess_vtk.jl")

end