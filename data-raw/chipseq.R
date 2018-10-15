# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next:
#     * Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package

library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(magrittr)
library(TeMPO)

ah <- AnnotationHub()

# Don't use whole genomes
# si <- SeqinfoForBSGenome("hg19")
# ss <- as(si, "GRanges")[c("chr20", "chr21", "chr22")]

# Only signal around CAGEclusters
data("CAGE_clusters")
ss <- reduce(CAGE_clusters + 400L)

query(ah, c("hela", "fc", "bigwig")) %>%
    mcols %>%
    View

H <- list(DNase="AH32877",
     H3K4me1="AH32879",
     H3K4me3="AH32881")#,
     #H3K27ac="AH32884")

gl <- lapply(H, function(x) ah[[x]]) %>%
    lapply(import.bw, which=ss) %>%
    lapply(subset, score > 0) %>%
    GRangesList()

bw_fnames <- file.path("inst", "extdata", paste0(names(gl), ".bw"))
mapply(export.bw, gl, bw_fnames)
