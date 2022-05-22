function volume(cell::Tetrahedron, nodes::Vector{Node{3,T}} where T <: Real, d::Vector{Float64})
    P = ones(Float64, 4, 4)
    for j = 1:4
        noded = getd(d, cell.nodes[j], 3)
        x = nodes[cell.nodes[j]].x
        for i = 1:3
            P[i+1,j] = x[i] + noded[i]
        end
    end
    return abs(det(P)) / 6
end

@inline function getd(d, node_id, dim)
    return d[(node_id-1)*dim+1:node_id*dim]
end