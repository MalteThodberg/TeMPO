context("Single or empty inputs")
library(TeMPO)
library(AnnotationHub)

#### SETUP ####

# Load data
data("CAGE_clusters")
data("CAGE_plus")

# Determine if bigwigs are available
ChIP_Seq <- list(DNase="AH32877",
                 H3K4me1="AH32879",
                 H3K4me3="AH32881",
                 H3K27ac="AH32884")

ah <- AnnotationHub()
fns <- na.omit(fileName(ah))

# Set flag accordingly
if(all(unlist(ChIP_Seq) %in% names(fns))){
    ChIP_Seq <- as(lapply(ChIP_Seq, function(x) ah[[x]]), "List")
    skip_bigwigs <- FALSE
}else{
    skip_bigwigs <- TRUE
}

#### TESTs ####

test_that("Single row on RleList input", {
    expect_s3_class(tidyMetaProfile(sites=CAGE_clusters[1],
                                    forward=CAGE_plus$WT),
                    class="tbl_df")
    expect_s3_class(tidyMetaProfile(sites=CAGE_clusters[1],
                                    forward=CAGE_plus$WT,
                                    reverse = CAGE_plus$KD),
                    class="tbl_df")
})

test_that("Single row on BigWig input", {
    # Skip if needed
    if(skip_bigwigs){skip("BigWigs not cached")}

    expect_s3_class(tidyMetaProfile(sites=CAGE_clusters[1],
                                    forward=ChIP_Seq$DNase),
                    class="tbl_df")
    expect_s3_class(tidyMetaProfile(sites=CAGE_clusters[1],
                                    forward=ChIP_Seq),
                    class="tbl_df")
})

test_that("Removing empty elements of GRangesList on RleList input", {
    expect_warning(tidyMetaProfile(sites=split(CAGE_clusters,
                                               CAGE_clusters$txType),
                                   forward=CAGE_plus$WT))
    expect_warning(tidyMetaProfile(sites=split(CAGE_clusters,
                                               CAGE_clusters$txType),
                                   forward=CAGE_plus$WT,
                                   reverse = CAGE_plus$KD))
})

test_that("Removing empty elements of GRangesList on BigWig input", {
    # Skip if needed
    if(skip_bigwigs){skip("BigWigs not cached")}

    # Gives warnings
    expect_warning(tidyMetaProfile(sites=split(CAGE_clusters,
                                               CAGE_clusters$txType),
                                   forward=ChIP_Seq$DNase))

    # BiocParallel does not pass on warning
    expect_s3_class(tidyMetaProfile(sites=split(CAGE_clusters,
                                               CAGE_clusters$txType),
                                   forward=ChIP_Seq),
                    class="tbl_df")
})
