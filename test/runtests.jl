using PhyloBuilder
using PhyloIO
using Random
using Test

const FIX = joinpath(@__DIR__, "fixtures")
fixture(name) = joinpath(FIX, name)
readmatrix(name) = permutedims(reduce(hcat, [parse.(Float64, split(l)) for l in eachline(fixture(name)) if !isempty(strip(l))]))

# patristic distances in a fixed taxon order: equal for two trees exactly when
# they have the same topology and branch lengths (all lengths > 0)
function pathdist(t, names)
    p = patristic(t)
    idx = [findfirst(==(n), p.names) for n in names]
    return p.d[idx, idx]
end

@testset "PhyloBuilder" begin
    include("alignment.jl")
    include("distances.jl")
    include("trees.jl")
end
