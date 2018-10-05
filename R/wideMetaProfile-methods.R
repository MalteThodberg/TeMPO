#' @rdname wideMetaProfile
#' @export
setMethod("wideMetaProfile",
          signature(forward="BigWigFile_OR_RleList", reverse="NULL_OR_missing",
                    sites="GenomicRanges"),
          function(sites, forward, reverse=NULL, upstream=100, downstream=100) {
              # Pre-checks
              assertthat::assert_that(all(width(sites) == 1),
                                      assertthat::is.count(upstream),
                                      assertthat::is.count(downstream))

              stopifnot(all(width(sites) == 1))

              # Expand
              sites <- promoters(sites, upstream=upstream,
                                 downstream=downstream)

              # Import
              o <- agnosticImport(signal=forward, sites=sites)

              # Reverse minus strand
              o <- revElements(o, strand(sites) == "-")

              # Coerce to matrix
              o <- do.call(rbind, o)

              # Add names
              colnames(o) <- c(seq(from=-upstream, to=-1),
                               seq_len(downstream))

              # Return
              o
          }
)

#' @rdname wideMetaProfile
#' @export
setMethod("wideMetaProfile",
          signature(forward="BigWigFile_OR_RleList",
                    reverse="BigWigFile_OR_RleList",
                    sites="GenomicRanges"),
          function(sites, forward, reverse=NULL, upstream=100, downstream=100) {
              # Pre-checks
              assertthat::assert_that(assertthat::are_equal(class(forward),
                                    class(reverse)),
                          all(width(sites) == 1),
                          assertthat::is.count(upstream),
                          assertthat::is.count(downstream))

              # Expand and split
              sites <- promoters(sites, upstream=upstream,
                                 downstream=downstream)
              strand(sites)[strand(sites) == "*"] <- "+"
              sites_stranded <- split(sites, strand(sites))

              # Coverage in all directions
              sense_forward <- agnosticImport(signal=forward,
                                              sites=sites_stranded$`+`)
              anti_forward <- agnosticImport(signal=reverse,
                                             sites=sites_stranded$`+`)
              sense_reverse <- agnosticImport(signal=reverse,
                                              sites=sites_stranded$`-`)
              anti_reverse <- agnosticImport(signal=forward,
                                             sites=sites_stranded$`-`)

              # Reverse minues
              sense_reverse <- revElements(x=sense_reverse,
                                           i=seq_along(sense_reverse))
              anti_reverse <- revElements(x=anti_reverse,
                                          i=seq_along(anti_reverse))

              # Combine matrices
              mat_sense <- do.call(rbind, c(sense_forward, sense_reverse))
              mat_anti <- do.call(rbind, c(anti_forward, anti_reverse))

              # Name
              colnames(mat_sense) <- colnames(mat_anti) <- c(seq(from=-upstream,
                                                                 to=-1),
                                                            seq_len(downstream))

              # Return
              list(sense=mat_sense, anti=mat_anti)
          })
