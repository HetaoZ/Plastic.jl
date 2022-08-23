# 1D and 3D
create_grid(c, nel::NTuple{dim,Int}, left::NTuple{dim,T}, right::NTuple{dim,T}) where dim where T = Ferrite.generate_grid(c, nel, Vec(left), Vec(right)) 

# 2D
create_grid(c, nel::NTuple{2,Int}, LL::NTuple{2,T}, LR::NTuple{2,T}, UR::NTuple{2,T}, UL::NTuple{2,T}) where T = Ferrite.generate_grid(c, nel, Vec(LL), Vec(LR), Vec(UR), Vec(UL))

function create_grid(c, nel::NTuple{2,Int}, LL::NTuple{2,T}, UR::NTuple{2,T}) where T 
    LR = (UR[1], LL[2])
    UL = (LL[1], UR[2])
    return create_grid(c, nel, LL, LR, UR, UL)
end