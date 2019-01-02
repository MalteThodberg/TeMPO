#' @rdname quantileTrim
#' @export
setMethod("quantileTrim",
          signature(mat1="matrix", mat2="NULL_OR_missing"),
          function(mat1, mat2=NULL, lower=0.0, upper=1.0) {
              # Pre-checks
              assertthat::assert_that(assertthat::not_empty(mat1))

              # Get rowsums
              x <- matrixStats::rowSums2(mat1)

              # Get limits
              lims <- quantile(x, c(lower, upper))

              # Find values within limits
              o <- x >= lims[1] & x <= lims[2]

              # Check if all sites were removed
              if(!any(o)){
                  stop("Excessive quantile trimming: ",
                       "All sites were removed!")
              }

              # Return filter
              mat1[o,,drop=FALSE]
          })

#' @rdname quantileTrim
#' @importFrom matrixStats rowSums2
#' @export
setMethod("quantileTrim",
          signature(mat1="matrix", mat2="matrix"),
          function(mat1, mat2=NULL, lower=0.0, upper=1.0) {
              # Pre-checks
              assertthat::assert_that(assertthat::not_empty(mat1),
                                      assertthat::not_empty(mat2),
                                      all(dim(mat1) == dim(mat2)))

              # Get rowsums
              x1 <- rowSums2(mat1)
              x2 <- rowSums2(mat2)

              # Get limits
              lims1 <- quantile(x1, c(lower, upper))
              lims2 <- quantile(x2, c(lower, upper))

              # Find values within limits
              o1 <- x1 >= lims1[1] & x1 <= lims1[2]
              o2 <- x2 >= lims2[1] & x2 <= lims2[2]
              o <- o1 & o2

              # Check if all sites were removed
              if(!any(o)){
                  stop("Excessive quantile trimming: ",
                       "All sites were removed!")
              }

              # Return filter
              list(mat1[o,,drop=FALSE],
                   mat2[o,,drop=FALSE])
          })
