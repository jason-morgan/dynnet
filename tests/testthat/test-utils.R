context("utility functions")

z <- matrix(c(1.0, 1.0,
              1.0, 0.0,
              0.5, 0.3), ncol=2, byrow=TRUE)

t <- matrix(c(1.0, 1.0,
              1.0, 0.0,
              0.5, 0.3), ncol=2, byrow=TRUE)

test_that("as_distance_vector returns correct structure", {
    expect_equal(class(as_distance_vector(z)), "matrix")
    expect_equal(dim(as_distance_vector(z))[2], 1)
    expect_equal(nrow(as_distance_vector(z)), nrow(z) * (nrow(z)-1) / 2)
})

test_that("as_distance_vector is correct for euclidean distances", {
    expect_equal(as_distance_vector(z)[1,1], 1.0)
})

test_that("calc_positions has proper structure", {
    expect_equal(length(calc_positions(1, z, t)), 1)
    expect_equal(length(calc_positions(1:2, z, t)), 2)
})

test_that("calc_positions returns correct results for euclidean distance", {
    expect_equal(calc_positions(1, z, t)[[1]][1,1], z[1,1])
    expect_equal(calc_positions(1:2, z, t)[[2]][1,1], z[1,1]*2)
})

lst1 <- list(c(1,2), c(1,2), c(1,2))
lst2 <- list(c(1,2), c(1,2), c(1))

test_that("list_all_equal returns the correct result", {
    expect_true(list_all_equal(lst1))
    expect_false(list_all_equal(lst2))
})

test_that("rksphere returns a vector of length 1", {
    expect_equal(sqrt(sum(rksphere(1)^2)), 1)
    expect_equal(sqrt(sum(rksphere(2)^2)), 1)
    expect_equal(sqrt(sum(rksphere(3)^2)), 1)
})

