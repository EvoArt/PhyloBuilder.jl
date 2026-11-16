# PhyloBuilder.jl

Phylogenetic trees from aligned nucleotide sequences. PhyloBuilder reads
alignments (FASTA, PHYLIP, Nexus), computes pairwise distances under the
usual substitution models (p, JC69, K80, F81, TN93, with optional gamma;
LogDet, paralinear) and builds UPGMA, WPGMA, neighbour-joining and BIONJ
trees. Distances and trees match R's `ape` and `phangorn`.

Trees are [`PhyloIO.PhyloTree`](https://github.com/EvoArt/PhyloIO.jl)s, so
they can be written as Newick or Nexus with PhyloIO and drawn with
[PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl).

## Installation

```julia
using Pkg
Pkg.add("PhyloBuilder")
```

## Quick start

```julia
using PhyloBuilder, PhyloIO

aln  = readalignment("aligned.fasta")
d    = distances(aln; model=:K80)
tree = nj(d; root=:midpoint)
writenewick("tree.nwk", tree)
```
