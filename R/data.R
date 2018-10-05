#' TeMPO example CAGE data
#'
#' Subset of publically available CAGE data from "Nuclear stability and
#' transcriptional directionality separate functionally distinct RNA species" by
#' Andersson et al.
#'
#' @format Three different data objects: \describe{
#'   \item{CAGE_clusters}{Location of TSSs and enhancers found with the
#'   CAGEfightR package} \item{CAGE_plus}{CAGE signal (Pooled Tags-Per-Million
#'   across each group) on the plus strand} \item{CAGE_minus}{CAGE signal
#'   (Pooled Tags-Per-Million across each group) on the minus strand} }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62047}
"CAGE_clusters"

#' @rdname CAGE_clusters
"CAGE_plus"

#' @rdname CAGE_clusters
"CAGE_minus"
