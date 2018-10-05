#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname agnosticImport
#' @export
setMethod("agnosticImport",
          signature(signal="BigWigFile", sites="GenomicRanges"),
          function(signal, sites) {
              # Import
              o <- suppressWarnings(import.bw(con=signal,
                                              which=sites,
                                              as="NumericList"))

              # Check all sites seqlevels are in signal seqlevels
              if(!any(seqlevels(signal) %in% seqlevels(sites))) {
                        stop("Seqlevel inconsistency between sites and signal")
              }

              # Check dimensions
              if(length(o) != length(sites)){
                  stop("BigWigFile was not properly imported!") # Output path?
              }

              # Pass on names
              names(o) <- names(sites)

              # Return
              o
          })

#' @rdname agnosticImport
#' @export
setMethod("agnosticImport",
          signature(signal="RleList", sites="GenomicRanges"),
          function(signal, sites) {
              # Subset and coerce
              o <- suppressWarnings(signal[sites])
              o <- methods::as(o, "NumericList")

              # Check all sites seqlevels are in signal seqlevels
              if(!any(seqlevels(signal) %in% seqlevels(sites))) {
                        stop("Seqlevel inconsistency between sites and signal")
              }

              # Check dimensions
              if(length(o) != length(sites)){
                  stop("RleList was not properly imported!")
              }

              # Pass on names
              names(o) <- names(sites)

              # Return
              o
          })
