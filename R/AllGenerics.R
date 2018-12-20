#' @include aaa.R
NULL

#### Helpers ####

#' Advanced usage: Import chunks of a genome-wide vector
#'
#' Import regions from coverage stored as a BigWigFile or RleList and return a NumericList.
#'
#' @param signal BigWigFile or RleList: Genomic signal
#' @param sites GenomicRanges: Regions to be extracted
#'
#' @return NumericList
#' @export
#'
#' @examples
#' data("CAGE_clusters")
#' data("CAGE_plus")
#'
#' # Import from an RleList
#' agnosticImport(signal=CAGE_plus$WT, sites=subset(CAGE_clusters, strand=="+"))
#'
#' # Import from a BigWig
#' \dontrun{
# # Get DNase data from AnnotationHub
# library(AnnotationHub)
# ah <- AnnotationHub()
# DNase <- ah[["AH32877"]]
#
# # Import only parts of the BigWig to conserve memory
# agnosticImport(signal=DNase, sites=subset(CAGE_clusters, strand=="+"))
#' }
setGeneric("agnosticImport", function(signal, sites) {
    standardGeneric("agnosticImport")
})

#' Advanced usage: Quantile trim a matrix or matrix pair
#'
#' Internal function for performing the quantile trimming. Removes outlier rows
#' with the highest or lowest sums. If two matrices are supplied, trims both
#' matrices simultaniously.
#'
#' @param mat1 matrix: Sites in rows, positions in columns.
#' @param mat2 matrix or NULL: Optional second matrix with the same shape.
#' @param lower numeric: Lower cutoff for trimming (between 0-1).
#' @param upper numeric: Upper cutoff for trimming (between 0-1).
#'
#' @return trimmed matrix or list of trimmed matrices
#' @export
#'
#' @examples
#' # Make some random matrices
#' m1 <- replicate(5, rnorm(10))
#' m2 <- replicate(5, rnorm(10))
#'
#' # Trim a single matrix
#' quantileTrim(mat1=m1, upper=0.75)
#'
#' # Trim a pair of matrices (e.g. forward/reverse strand)
#' quantileTrim(mat1=m1, mat2=m2, upper=0.75)
setGeneric("quantileTrim",
           function(mat1, mat2=NULL, lower=0.0, upper=1.0){
               standardGeneric("quantileTrim")
           })

#### Main functions ####

#' Advanced usage: Meta Profile in wide format
#'
#' Calculates the wide-format meta profile, with regions in rows and positions in columns.
#'
#' @param sites GenomicRanges: Single-bp genomic locations to calculate the
#'   meta-profile over.
#' @param forward BigWigFile or RleList: Genome-wide signal stored either on
#'   disk as a BigWig-file or in memory as an RleList-object.
#' @param reverse BigWigFile, RleList or NULL: If reverse=NULL, the forward
#'   signal is taken as unstranded. If not the genomic signal is taken to be
#'   stranded. The class of reverse must be the same as forward.
#' @param upstream integer: Number of bases to extend upstream.
#' @param downstream integer: Number of bases to extend downstream
#'
#' @return if reverse=NULL a matrix, else a list of two matrices with names
#'   "sense" and "anti"
#' @export
#'
#' @examples
#' data("CAGE_clusters")
#' data("CAGE_plus")
#' data("CAGE_minus")
#'
#' enhancers <- subset(CAGE_clusters, clusterType == "enhancer")
#'
#' # Unstranded data gives a single matrix:
#' wideMetaProfile(sites=enhancers,
#'                 forward=CAGE_plus$WT)
#'
#' # Stranded data gives two matrices (sense and anti):
#' wideMetaProfile(sites=enhancers,
#'                 forward=CAGE_plus$WT,
#'                 reverse=CAGE_minus$WT)
setGeneric("wideMetaProfile", function(sites, forward, reverse=NULL,
                                       upstream=100, downstream=100) {
    standardGeneric("wideMetaProfile")
    })

#' TeMPO: Tidy Meta Profiles
#'
#' @inheritParams wideMetaProfile
#' @param trimLower numeric: Lower quantile used for trimming.
#' @param trimUpper numeric: Upper quantile used for trimming.
#' @param sumFun function: The function used for summarising the metaprofile.
#'   Should be similar to colSums.
#'
#' @return A tibble with the following columns
#' \describe{
#'   \item{pos0}{0-based position relative to sites, useful for plotting}
#'   \item{pos1}{1-based position relative to sites}
#'   \item{sense}{Summarized signal for each position}
#' }
#' Depending on the input, addditional columns might be added
#' \describe{
#'   \item{anti}{If reverse is supplied the genome-wide signal is assumed to be stranded, and the output is split into sense/antisense relative to the orientation of sites.}
#'   \item{sites}{If sites is a GRangesList, indicates the different sets of sites}
#'   \item{signal}{if forward/reverse is a BigWigFileList or an RleListList, indicates the different genomic signals}
#' }
#' @export
#'
#' @examples
#' data("CAGE_clusters")
#' data("CAGE_plus")
#' data("CAGE_minus")
#'
#' # Stranded Meta-profile across enhancers:
#' P <- tidyMetaProfile(sites = subset(CAGE_clusters, clusterType=="enhancer"),
#'                      forward = CAGE_plus,
#'                      reverse = CAGE_minus,
#'                      upstream=200,
#'                      downstream=200)
#'
#' # Plot the resulting tibble with ggplot2:
#' library(tidyr)
#' library(ggplot2)
#'
#' P %>%
#' gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
#' ggplot(aes(x=pos0, y=score, color=direction)) +
#' geom_line(alpha=0.75) +
#' scale_color_brewer("Direction", palette="Set1") +
#' labs(x="Relative position from center",
#' y="Average Signal")
#'
#' # See the TeMPO vignette for more examples and settings!
setGeneric("tidyMetaProfile",
           function(sites, forward, reverse=NULL,
                    upstream=100, downstream=100,
                    trimLower=0.0, trimUpper=1.0,
                    sumFun=matrixStats::colMeans2){
               standardGeneric("tidyMetaProfile")
               })
