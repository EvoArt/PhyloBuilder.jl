"""
    PhyloBuilder

Build phylogenies from aligned nucleotide sequences: distances, UPGMA,
WPGMA, neighbour joining and BIONJ (more to come: parsimony, maximum
likelihood, bootstrap support). Trees are `PhyloIO.PhyloTree`s, so they
can be written with PhyloIO and drawn with PlotPhyloTree.
"""
module PhyloBuilder

using Printf: Printf, @sprintf
using Random: Random, AbstractRNG
using PhyloIO
using PhyloIO: PhyloTree, addnode!, nnodes, ntips, isleaf, isroot, parentnode, childnodes, branchlength,
               nodelabel, rootnode, preorder, tips, mrca, nodeindex, rootdistance, ladderize!, reroot

export Alignment, ntaxa, nsites, taxa, sequences, basefreqs,
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
