#' @rdname tidyMetaProfile
#' @importFrom tibble tibble
#' @export
setMethod("tidyMetaProfile",
          signature(forward="BigWigFile_OR_RleList",
                    reverse="NULL_OR_missing",
                    sites="GenomicRanges"),
          function(sites, forward, reverse=NULL, upstream=100, downstream=100,
                   trimLower=0, trimUpper=1, sumFun=colMeans) {
              # Pre-checks
              assertthat::assert_that(assertthat::is.scalar(trimLower),
                                      assertthat::is.scalar(trimUpper),
                                      trimLower >= 0 & trimLower < 1,
                                      trimUpper > 0 & trimUpper <= 1,
                                      is.function(sumFun))

              # Import
              o <- wideMetaProfile(sites=sites,
                                   forward=forward,
                                   reverse=NULL,
                                   upstream=upstream,
                                   downstream=downstream)

              # Trim
              o <- quantileTrim(o, upper=trimUpper, lower=trimLower)

              # Summarise
              o <- sumFun(o)

              # Tidy up
              o <- tibble(pos1=c(seq(from=-upstream, to=-1),
                                         seq_len(downstream)),
                                  pos0=seq(from=-upstream, to=downstream-1),
                                  sense=o)

              # Return
              o
          })

#' @rdname tidyMetaProfile
#' @importFrom tibble tibble
#' @export
setMethod("tidyMetaProfile",
          signature(forward="BigWigFile_OR_RleList",
                    reverse="BigWigFile_OR_RleList",
                    sites="GenomicRanges"),
          function(sites, forward, reverse=NULL, upstream=100, downstream=100,
                   trimLower=0, trimUpper=1, sumFun=matrixStats::colMeans2) {
              # Pre-checks
              assertthat::assert_that(assertthat::is.scalar(trimLower),
                                      assertthat::is.scalar(trimUpper),
                                      trimLower >= 0 & trimLower < 1,
                                      trimUpper > 0 & trimUpper <= 1,
                                      is.function(sumFun))


              # Import
              o <- wideMetaProfile(sites=sites,
                                   forward=forward,
                                   reverse=reverse,
                                   upstream=upstream,
                                   downstream=downstream)

              # Trim
              o <- quantileTrim(mat1=o$sense, mat2=o$anti,
                                upper=trimUpper, lower=trimLower)

              # Summarise
              o_sense <- sumFun(o[[1]])
              o_anti <- sumFun(o[[2]])

              # Tidy up
              o <- tibble(pos1=c(seq(from=-upstream, to=-1),
                                         seq_len(downstream)),
                                  pos0=seq(from=-upstream, to=downstream-1),
                                  sense=o_sense,
                                  anti=o_anti)

              # Return
              o
          })

#' @rdname tidyMetaProfile
#' @importFrom dplyr bind_rows
#' @export
setMethod("tidyMetaProfile",
          signature(forward="BigWigFile_OR_RleList",
                    reverse="NULL_OR_missing",
                    sites="GenomicRangesList"),
          function(sites, forward, reverse=NULL, upstream=100, downstream=100,
                   trimLower=0, trimUpper=1, sumFun=matrixStats::colMeans2) {
              # Loop over elements
              o <- lapply(sites, tidyMetaProfile,
                          forward=forward, reverse=reverse,
                          upstream=upstream, downstream=downstream,
                          trimLower=trimLower, trimUpper=trimUpper,
                          sumFun=sumFun)

              # Bind rows
              o <- bind_rows(o, .id="sites")

              # Return
              o
          })

#' @rdname tidyMetaProfile
#' @importFrom dplyr bind_rows
#' @export
setMethod("tidyMetaProfile",
          signature(forward="BigWigFile_OR_RleList",
                    reverse="BigWigFile_OR_RleList",
                    sites="GenomicRangesList"),
          function(sites, forward, reverse=NULL,
                   upstream=100, downstream=100,
                   trimLower=0, trimUpper=1, sumFun=matrixStats::colMeans2) {
              # Loop over elements
              o <- lapply(sites, tidyMetaProfile,
                          forward=forward, reverse=reverse,
                          upstream=upstream, downstream=downstream,
                          trimLower=trimLower, trimUpper=trimUpper,
                          sumFun=sumFun)

              # Bind rows
              o <- bind_rows(o, .id="sites")

              # Return
              o
          })

#' @rdname tidyMetaProfile
#' @importFrom dplyr bind_rows
#' @export
setMethod("tidyMetaProfile",
          signature(forward="BigWigFileList_OR_RleListList",
                    reverse="NULL_OR_missing",
                    sites="GenomicRanges_OR_GenomicRangesList"),
          function(sites, forward, reverse=NULL,
                   upstream=100, downstream=100,
                   trimLower=0, trimUpper=1, sumFun=matrixStats::colMeans2) {
              # Loop over elements
              o <- bplapply(forward, tidyMetaProfile,
                            sites=sites, reverse=reverse,
                            upstream=upstream, downstream=downstream,
                            trimLower=trimLower, trimUpper=trimUpper,
                            sumFun=sumFun)

              # Bind rows
              o <- bind_rows(o, .id="signal")

              # Return
              o
          })

#' @rdname tidyMetaProfile
#' @importFrom dplyr bind_rows
#' @export
setMethod("tidyMetaProfile",
          signature(forward="BigWigFileList_OR_RleListList",
                    reverse="BigWigFileList_OR_RleListList",
                    sites="GenomicRanges_OR_GenomicRangesList"),
          function(sites, forward, reverse=NULL,
                   upstream=100, downstream=100,
                   trimLower=0, trimUpper=1, sumFun=matrixStats::colMeans2) {
              # Loop over elements
              o <- bpmapply(tidyMetaProfile, forward, reverse,
                          MoreArgs=list(sites=sites,
                                        upstream=upstream,
                                        downstream=downstream,
                                        trimLower=trimLower,
                                        trimUpper=trimUpper,
                                        sumFun=sumFun),
                          SIMPLIFY = FALSE)

              # Bind rows
              o <- bind_rows(o, .id="signal")

              # Return
              o
          })
