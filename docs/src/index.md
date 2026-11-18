# PhyloBuilder.jl

Phylogenetic trees from aligned nucleotide sequences. PhyloBuilder reads
alignments (FASTA, PHYLIP, Nexus), computes pairwise distances under the
usual substitution models (p, JC69, K80, F81, TN93, with optional gamma;
LogDet, paralinear) and builds UPGMA, WPGMA, neighbour-joining and BIONJ
trees. Distances and trees match R's `ape` and `phangorn`.

Trees use a temporary bundled copy of [PhylogenyIO](https://github.com/EvoArt/PhylogenyIO.jl).
Use `PhyloBuilder.PhylogenyIO` for Newick and Nexus I/O. Its `PhyloTree` type
differs from the copy bundled by [PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl),
so transfer trees through a Newick string or file until both packages share the
registered library.

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
