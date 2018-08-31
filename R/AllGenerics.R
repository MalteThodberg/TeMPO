#### Helpers ####

#' Advanced usage: Import chunks of a genome-wide vector
#'
#' @param signal BigWigFile or RleList
#' @param sites GenomicRanges
#'
#' @return NumericList
#' @export
#'
#' @examples
setGeneric("agnosticImport", function(signal, sites) {
    standardGeneric("agnosticImport")
})

#' Advanced usage: Quantile trim a matrix or matrix pair
#'
#' @param mat1 matrix
#' @param mat2 matrix or NULL
#' @param lower numeric
#' @param upper numeric
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
#'
#' @param sites GenomicRanges: Single-bp genomic locations to calculate the meta-profile over.
#' @param forward BigWigFile or RleList: Genome-wide signal stored either on disk as a BigWig-file or in memory as an RleList-object.
#' @param reverse BigWigFile, RleList or NULL: If reverse=NULL, the forward signal is taken as unstranded. If not the genomic signal is taken to be stranded. The class of reverse must be the same as forward.
#' @param upstream integer: Number of bases to extend upstream.
#' @param downstream integer: Number of bases to extend downstream
#'
#' @return if reverse=NULL a matrix, else a list of two matrices with names "sense" and "anti"
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
#' @param sumFun function: The function used for summarising the metaprofile. Should be similar to colSums.
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
setGeneric("tidyMetaProfile",
           function(sites, forward, reverse=NULL,
                    upstream=100, downstream=100,
                    trimLower=0.0, trimUpper=1.0,
                    sumFun=matrixStats::colMeans2){
               standardGeneric("tidyMetaProfile")
               })
