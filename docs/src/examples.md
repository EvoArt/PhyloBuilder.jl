# Examples

## Distances and trees

A small alignment of eight sequences shipped with the tests:

```@example main
using PhyloBuilder, PhyloIO

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

p = treeplot(t_nj)
scalebar!(p)
p
```

```julia
savefig("tree.svg", p)
```

## From ASVs to a tree

Amplicon sequence variants from
[DenoiseDNA](https://github.com/EvoArt/DenoiseDNA.jl), aligned with MAFFT
(via `MAFFT_jll`), then a BIONJ tree drawn with PlotPhyloTree. `seqtab_nochim`
is the chimera-free sequence table from the DenoiseDNA tutorial.

!!! note
    `MAFFT_jll` has no Windows build, so this example runs on Linux and macOS
    only. A pure-Julia aligner, MSABuilder.jl, is in development.

```julia
using DenoiseDNA, PhyloBuilder, PlotPhyloTree, MAFFT_jll

asvs = seqtab_nochim.sequences
open("asvs.fasta", "w") do io
    for (i, s) in enumerate(asvs)
        println(io, ">ASV", i, '\n', s)
    end
end
run(pipeline(`$(mafft()) --auto --quiet asvs.fasta`; stdout="asvs_aligned.fasta"))

aln  = readalignment("asvs_aligned.fasta")
tree = bionj(aln; model=:TN93, gamma=0.5, root=:midpoint, ladderize=true)

p = treeplot(tree; layout=:circular)
savefig("asv-tree.svg", p)
```
