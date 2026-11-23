# PhyloBuilder.jl

Phylogenetic trees from aligned nucleotide sequences. PhyloBuilder reads
alignments (FASTA, PHYLIP, Nexus), computes pairwise distances under the
usual substitution models (p, JC69, K80, F81, TN93, with optional gamma;
LogDet, paralinear) and builds UPGMA, WPGMA, neighbour-joining and BIONJ
trees. Distances and trees match R's `ape` and `phangorn`.

Trees are [PhylogenyIO](https://github.com/EvoArt/PhylogenyIO.jl) `PhyloTree`s.
`PhyloBuilder.PhylogenyIO` gives Newick and Nexus I/O without a separate install,
and [PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl) draws the trees
directly.

## Installation

```julia
using Pkg
Pkg.add("PhyloBuilder")
```

## Quick start

```julia
using PhyloBuilder, PhyloBuilder.PhylogenyIO

aln  = readalignment("aligned.fasta")
d    = distances(aln; model=:K80)
tree = nj(d; root=:midpoint)
writenewick("tree.nwk", tree)
```
