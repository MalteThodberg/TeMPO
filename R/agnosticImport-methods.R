#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname agnosticImport
#' @export
setMethod("agnosticImport",
          signature(signal = "BigWigFile", sites = "GenomicRanges"),
          function(signal, sites) {
              # Pre-checks
              # Necessary?

              # Remove Out-of-Bounds
              si <- seqinfo(signal)
              sites <- remove_OoB(sites, si)

              # Import
              o <- import.bw(con = signal,
                             which = sites,
                             as = "NumericList")

              # Check dimensions - add path here?
              if (length(o) != length(sites)) {
                  stop("Unknown error: ",
                       "Not all sites were imported from the BigWig") # Path?
              }

              if (dplyr::n_distinct(elementNROWS(o)) > 1) {
                  stop("Unknown error: ",
                       "Sites were not correctly imported from the BigWig")
              }

              # Pass on names
              names(o) <- names(sites)

              # Return
              o
          })

#' @rdname agnosticImport
#' @export
setMethod("agnosticImport",
          signature(signal = "RleList", sites = "GenomicRanges"),
          function(signal, sites) {
              # Pre-checks
              # Necessary?

              # Remove Out-of-Bounds
              si <- elementNROWS(signal)
              si <- Seqinfo(seqnames = names(si),
                            seqlengths = si)
              sites <- remove_OoB(sites, si)

              # Import
              o <- signal[sites]
              o <- methods::as(o, "NumericList")

              # Check dimensions
              if (length(o) != length(sites)) {
                  stop("Unknown error: ",
                       "Not all sites were imported from the RleList") # Path?
              }

              if (dplyr::n_distinct(elementNROWS(o)) > 1) {
                  stop("Unknown error: ",
                       "Sites were not correctly imported from the RleList")
              }

              # Pass on names
              names(o) <- names(sites)

              # Return
              o
          })
