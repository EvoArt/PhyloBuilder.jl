# PhyloBuilder.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/PhyloBuilder.jl/dev)
[![CI](https://github.com/EvoArt/PhyloBuilder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/PhyloBuilder.jl/actions/workflows/CI.yml)

Phylogenetic trees from aligned nucleotide sequences: pairwise distances
(p, JC69, K80, F81, TN93 with optional gamma; LogDet, paralinear) and
UPGMA, WPGMA, neighbour-joining and BIONJ trees. Results match R's `ape` and
`phangorn`. This release bundles the `PhylogenyIO` tree implementation.
`PhyloBuilder.PhylogenyIO.PhyloTree` is distinct from the tree type bundled by
[PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl). Transfer trees through
Newick until both packages depend on the registered `PhylogenyIO` release.

```julia
using Pkg
Pkg.add("PhyloBuilder")
```

```julia
using PhyloBuilder, PhyloBuilder.PhylogenyIO

aln  = readalignment("aligned.fasta")        # FASTA, PHYLIP or Nexus
d    = distances(aln; model=:TN93, gamma=0.5)
tree = bionj(d; root=:midpoint)
writenewick("tree.nwk", tree)
```

## Benchmarks

500 sequences × 1500 sites, times in ms (Intel i7-8665U).

| step | PhyloBuilder (1 thread) | PhyloBuilder (8 threads) | R |
|---|---:|---:|---:|
| distances, JC69 | 38 | 21 | 470 (`ape::dist.dna`) |
| distances, K80 | 40 | 21 | 1110 |
| distances, TN93 | 43 | 21 | 1300 |
| distances, LogDet | 42 | 22 | 2000 |
| NJ | 61 | – | 95 (`ape::nj`) |
| BIONJ | 60 | – | 235 (`ape::bionj`) |
| UPGMA | 6 | – | 7 (`phangorn::upgma`) |
