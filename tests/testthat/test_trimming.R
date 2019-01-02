context("Excessive trimming")
m1 <- replicate(5, rnorm(10))
m2 <- replicate(5, rnorm(10))
m3 <- matrix(ncol = 5)

test_that("Empty or misshapen input", {
    expect_error(quantileTrim(m3, lower = 0, upper = 1))
    expect_error(quantileTrim(m1, m3, lower = 0, upper = 1))

})

test_that("Everything trimmed away", {
    expect_error(quantileTrim(m1, lower = 1, upper = 0))
    expect_error(quantileTrim(m1, m2, lower = 1, upper = 0))
})
