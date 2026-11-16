# Reference values from ape/phangorn for the PhyloBuilder tests.
# Run from the package root:  Rscript test/fixtures/make_fixtures.R
# small.fa: 8 taxa × 300 sites with a gap run and R/N codes.
suppressMessages({library(ape); library(phangorn)})
dir <- file.path("test", "fixtures")
x <- read.dna(file.path(dir, "small.fa"), format = "fasta")

wmat <- function(m, name)
    write.table(formatC(unname(as.matrix(m)), digits = 17, format = "g"),
                file.path(dir, name), quote = FALSE, row.names = FALSE, col.names = FALSE)

cases <- list(raw = list("N"), p = list("raw"), JC69 = list("JC69"), K80 = list("K80"),
              F81 = list("F81"), TN93 = list("TN93"), logdet = list("logdet"),
              paralin = list("paralin"),
              JC69_g05 = list("JC69", 0.5), K80_g05 = list("K80", 0.5),
              F81_g07 = list("F81", 0.7), TN93_g07 = list("TN93", 0.7))
for (nm in names(cases)) {
    a <- cases[[nm]]
    d <- if (length(a) == 1) dist.dna(x, model = a[[1]], pairwise.deletion = TRUE, as.matrix = TRUE)
         else dist.dna(x, model = a[[1]], gamma = a[[2]], pairwise.deletion = TRUE, as.matrix = TRUE)
    wmat(d, paste0("dist_", nm, ".txt"))
}
writeLines(sprintf("%.17g", base.freq(x)), file.path(dir, "basefreq.txt"))

d <- dist.dna(x, model = "K80", pairwise.deletion = TRUE)
write.tree(nj(d), file.path(dir, "nj_K80.nwk"), digits = 17)
write.tree(bionj(d), file.path(dir, "bionj_K80.nwk"), digits = 17)   # ape's BIONJ works in single precision
write.tree(phangorn::upgma(d), file.path(dir, "upgma_K80.nwk"), digits = 17)
write.tree(phangorn::wpgma(d), file.path(dir, "wpgma_K80.nwk"), digits = 17)
