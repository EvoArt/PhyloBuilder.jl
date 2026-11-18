# Examples

## Distances and trees

A small alignment of eight sequences shipped with the tests:

```@example main
using PhyloBuilder, PhyloBuilder.PhylogenyIO

aln = readalignment(pkgdir(PhyloBuilder, "test", "fixtures", "small.fa"))
ntaxa(aln), nsites(aln)
```

Pairwise distances under the Kimura two-parameter model, with gamma-distributed
rates:

```@example main
d = distances(aln; model=:K80, gamma=0.5)
round.(d; digits=3)
```

```@example main
d["s1", "s2"]
```

Trees from the distance matrix, or straight from the alignment:

```@example main
t_nj    = nj(d; root=:midpoint, ladderize=true)
t_bionj = bionj(aln; model=:TN93, root="s1")
t_upgma = upgma(d)
writenewick(t_nj)
```

## Drawing a tree

With [PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl):

```@example main
using PlotPhyloTree

p = treeplot(writenewick(t_nj))
scalebar!(p)
p
```

```julia
savefig("tree.svg", p)
```

## From ASVs to a tree

Amplicon sequence variants from
[DenoiseDNA](https://github.com/EvoArt/DenoiseDNA.jl), aligned with
[MSABuilder](https://github.com/EvoArt/MSABuilder.jl), then a tree drawn with
PlotPhyloTree. `seqtab_nochim` is the chimera-free sequence table from the
DenoiseDNA tutorial.

```julia
using DenoiseDNA, MSABuilder, PhyloBuilder, PlotPhyloTree

tree = seqtab_nochim |> align |> nj        # ASVs named ASV1, ASV2, …
p = treeplot(midpoint(tree); layout=:circular)
savefig("asv-tree.svg", p)
```

Any collection of `name => sequence` pairs can be passed where an alignment
is expected, so `nj(align(seqs))` works for plain sequences too.
