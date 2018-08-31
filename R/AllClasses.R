#### New classes ####

#' @export
#' @rdname  RleListList
.RleListList <- setClass("RleListList",
                         contains="SimpleList",
                         prototype=prototype(elementType="RleList"))

#' List of RleList-objects
#'
#' @param ... RleLists
#'
#' @return RleListList
#' @export
#'
#' @examples
RleListList <- function(...){
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
              members = c("BigWigFile", "RleList"))) # Gives warning
setClassUnion("NULL_OR_missing",
              members = c("NULL", "missing"))
setClassUnion("BigWigFileList_OR_RleListList",
              members = c("BigWigFileList", "RleListList"))
