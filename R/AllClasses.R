#' @include aaa.R
NULL

#### New classes ####

#' @export
#' @rdname  RleListList
setClass("RleListList",
         contains = "SimpleList",
         prototype = prototype(elementType = "RleList"))

#' List of RleList-objects
#'
#' The RleListList is simply a List where all elements are RleList. The most
#' likely use case is a collection of genome-wide signals calculated using e.g.
#' coverage() or imported using import.bw().
#'
#' @param ... RleList-objects to be gathered into a list.
#'
#' @return RleListList
#' @export
#' @examples
#' # Example coverage calculate
#' gr <- GRanges(seqnames = paste0("chr1", 1:10),
#'               ranges = IRanges(1:10, width = 1))
#' cvg <- coverage(gr)
#'
#' # Use constructor:
#' RleListList(cvg, cvg)
#'
#' # Or simply coerce to a List
#' as(list(cvg, cvg), "List")
RleListList <- function(...) {
    # Simply coercion
    o <- List(...)

    # Check if the result is valid
    stopifnot(methods::is(o, "RleListList"))
    methods::validObject(o)

    # Return
    o
}

#### Class Unions ####

suppressWarnings(setClassUnion("BigWigFile_OR_RleList",
                               members = c("BigWigFile",
                                           "RleList"))) # Gives warning
setClassUnion("NULL_OR_missing",
              members = c("NULL", "missing"))
setClassUnion("BigWigFileList_OR_RleListList",
              members = c("BigWigFileList", "RleListList"))
