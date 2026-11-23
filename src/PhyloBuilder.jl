"""
    PhyloBuilder

Phylogenetic trees from aligned nucleotide sequences: pairwise distances,
UPGMA, WPGMA, neighbour joining and BIONJ. Trees are `PhylogenyIO.PhyloTree`s.
"""
module PhyloBuilder

using Printf: Printf, @sprintf
using Random: Random, AbstractRNG
using PhylogenyIO
using PhylogenyIO: PhyloTree, addnode!, nnodes, ntips, isleaf, isroot, parentnode, childnodes, branchlength,
               nodelabel, rootnode, preorder, tips, mrca, nodeindex, rootdistance, ladderize!, reroot

export PhylogenyIO

# `sequences` is not exported (CommunityEcology exports a different one); use
# `PhyloBuilder.sequences(aln)` or `using PhyloBuilder: sequences`.
export Alignment, ntaxa, nsites, taxa, basefreqs,
       readalignment, readfasta, readphylip, writefasta, writephylip,
       DistanceMatrix, distances, patristic,
       upgma, wpgma, nj, bionj,
       midpoint, unroot, root

include("style.jl")
include("progress.jl")
include("alignment.jl")
include("distances.jl")
include("trees.jl")
include("distancetrees.jl")

end
