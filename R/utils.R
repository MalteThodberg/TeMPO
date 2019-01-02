remove_empty <- function(grl){
    # Check
    stopifnot(methods::is(grl, "GRangesList"))

    # Find empty elements
    elements <- elementNROWS(grl) != 0

    # Remove if need
    if(!all(elements)){
        warning("Removed empty elements from GRangesList: ",
                sum(!elements), call. = FALSE)
        grl <- grl[elements]
    }

    # Return
    grl
}

remove_OoB <- function(gr, si){
    # Checks
    stopifnot(methods::is(gr, "GenomicRanges"),
              methods::is(si, "Seqinfo"))

    # Coerce
    si <- methods::as(si, "GRanges")

    # Overlap
    o <- suppressWarnings(subsetByOverlaps(x=gr,
                                           ranges=si,
                                           type="within",
                                           ignore.strand=TRUE))

    # Report
    l0 <- length(gr)
    l <- length(o)

    if(l == 0 & l0 != 0){
        stop("No sites could be imported: ",
             "Check the genome information of input objects with seqinfo()")
    }else if(l != l0){
        warning(l0 - l,
                " sites were removed on import due to incompatible genomes!",
                "Check the genome information of input objects with seqinfo()")
    }else if(l0 > l){
        stop("Too many sites were imported: ",
             "Check the genome information of input objects with seqinfo()")
    }

    # Clean seqlevels
    seqlevels(o) <- seqlevelsInUse(o)

    # Return
    o
}
