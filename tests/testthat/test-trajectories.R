context("trajectory functions")

z <- matrix(c(1.0, 1.0,
              1.0, 0.0,
              0.5, 0.3), ncol=2, byrow=TRUE)

t <- matrix(c(1.0, 1.0,
              1.0, 0.0,
              0.5, 0.3), ncol=2, byrow=TRUE)

test_that("trajectory_linear returns correct results", {
    expect_equal(trajectory_linear(1, z, t), z)
    expect_equal(trajectory_linear(2, z, t), z+t)    
})
