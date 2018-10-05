#' @rdname quantileTrim
#' @export
setMethod("quantileTrim",
          signature(mat1="matrix", mat2="NULL_OR_missing"),
          function(mat1, mat2=NULL, lower=0.0, upper=1.0) {
              # Get rowsums
              x <- matrixStats::rowSums2(mat1)

              # Get limits
              lims <- quantile(x, c(lower, upper))

              # Find values within limits
              o <- x >= lims[1] & x <= lims[2]

              # Return filter
              mat1[o,]
          })

#' @rdname quantileTrim
#' @importFrom matrixStats rowSums2
#' @export
setMethod("quantileTrim",
          signature(mat1="matrix", mat2="matrix"),
          function(mat1, mat2=NULL, lower=0.0, upper=1.0) {
              stopifnot(all(dim(mat1) == dim(mat2)))

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

              # Return filter
              list(mat1[o,], mat2[o,])
          })
