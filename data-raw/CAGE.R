library(CAGEfightR)
library(purrr)
library(magrittr)
library(TeMPO)

# Cnvert BE
if(FALSE){
    # Make file names
    bed_fnames <- list.files("~/Documents/Work/RStudio/TeMPO_data/GSE62047_RAW", full.names = TRUE)

    plus_fnames <- bed_fnames %>%
        stringr::str_replace(".bed.gz", "_plus.bw") %>%
        stringr::str_replace("GSE62047_RAW", "bw")

    minus_fnames <- bed_fnames %>%
        stringr::str_replace(".bed.gz", "_minus.bw") %>%
        stringr::str_replace("GSE62047_RAW", "bw")

    # Import and strand
    bed_grs <- map(bed_fnames, import, genome=SeqinfoForBSGenome("hg19"))
    bed_stranded <- map(bed_grs, CAGEfightR:::splitByStrand)
    bed_plus <- map(bed_stranded, ~.x$`+`)
    bed_minus <- map(bed_stranded, ~.x$`-`)

    # Export
    walk2(bed_plus, plus_fnames, export.bw)
    walk2(bed_minus, minus_fnames, export.bw)
}

#### Pooled ####

bw_files <- list.files("~/Documents/Work/RStudio/TeMPO_data/bw", full.names = TRUE)
bw_plus <- stringr::str_subset(bw_files, "_plus.bw") %>% BigWigFileList()
bw_minus <- stringr::str_subset(bw_files, "_minus.bw") %>% BigWigFileList()
names(bw_plus) <- names(bw_minus) <- c(paste0("Egfp.R", 1:3), paste0("Rrp40.R", 1:3))

register(SerialParam())

CTSSs <- quantifyCTSSs(bw_plus, bw_minus, genome=SeqinfoForBSGenome("hg19"))

WT <- CTSSs[,1:3] %>%
    calcTPM %>%
    calcPooled
KD <- CTSSs[,4:6] %>%
    calcTPM %>%
    calcPooled

#### Clusters ####

TCs <- quickTSSs(WT)
BCs <- quickEnhancers(WT)

#
TSSs <- TCs %>%
    calcTPM %>%
    subsetBySupport(inputAssay="TPM", unexpressed=1, minSamples=2) %>%
    rowRanges
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

TSSs <- assignTxType(TSSs, txdb, swap="thick")
# TSSs <- calcShape(TSSs, pooled=WT, outputColumn='IQR',
#                   shapeFunction=shapeIQR, lower=0.2, upper=0.8)

Enhancers <- rowRanges(BCs)
Enhancers <- assignTxType(Enhancers, txdb)
Enhancers <- subset(Enhancers, txType %in% c("intergenic", "intron"))

# Merge
Enhancers$clusterType <- "enhancer"
Enhancers$balance <- NULL
Enhancers$bidirectionality <- NULL
TSSs$clusterType <- "TSS"
TSSs$support <- NULL

CAGE_clusters <- c(Enhancers, TSSs[!overlapsAny(TSSs, Enhancers)]) %>%
    subset(seqnames %in% paste0("chr", 20:22)) %>%
    swapRanges() %>%
    subset(select=c(score, clusterType, txType))

#### Subset coverage ####

FS <- . %>%
    rowRanges %>%
    as("GRanges") %>%
    subsetByOverlaps(CAGE_clusters + 500) %>%
    #subset(seqnames %in% paste0("chr", 20:22)) %>%
    CAGEfightR:::splitByStrand() %>%
    lapply(coverage, weight="score") %>%
    lapply(round, digits=8)

WT_stranded <- WT %>% FS
KD_stranded <- KD %>% FS

CAGE_plus <- List(WT=WT_stranded$`+`, KD=KD_stranded$`+`)
CAGE_minus <- List(WT=WT_stranded$`-`, KD=KD_stranded$`-`)

#### Annotated ####
#
# # Extract promoters
# anno <- transcripts(txdb)
# anno$thick <- ranges(resize(anno, width = 1, fix="start"))
# anno <- promoters(anno, 50, 50)
#
# # Remove overlapping
# anno <- subsetByOverlaps(anno, disjoin(anno) %>% subset(width==100))
#
# # Quantify
# anno <- quantifyClusters(WT, anno)
#
# anno <- txdb %>%
#     promoters(50, 50) %>%
#     GenomicRanges::reduce() %>%
#     subset(seqnames %in% paste0("chr", 20:22) & width == 100)

#### Write to files ####

devtools::use_data(CAGE_clusters, overwrite = TRUE)
devtools::use_data(CAGE_plus, overwrite = TRUE)
devtools::use_data(CAGE_minus, overwrite = TRUE)

#
#
#
# combineClusters(object1 = TSSs, object2 = Enhancers, removeIfOverlapping="object1")
#
# # Types and shapes
# TSSs <- assignTxType(TSSs, TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#
# combineClusters(object1 = TSSs, object2 = Enhancers, removeIfOverlapping="object1")
#
# # Write file
#
