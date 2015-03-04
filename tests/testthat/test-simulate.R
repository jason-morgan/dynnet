context("network simulation")

## Number of nodes in each group.
n <- 20
n_1 <- split_n(n, 2)
n_2 <- split_n(n, 3)

## Positions
p_1 <- matrix(c(-2, 2))
p_2 <- matrix(c(-2, 0, 2, -2, 2, -2), ncol=2)

## Sigmas
s_1 <- list(diag(1), diag(1))
s_2 <- list(diag(2), diag(2), diag(2))

test_that("generate_group produces the correct structure", {
    expect_equal(dim(generate_group(10, 0, diag(1))), c(10, 1))
    expect_equal(dim(generate_group(10, c(0, 0), diag(2))), c(10, 2))
})

test_that("generate_latent_values produces the correct structure", {
    expect_equal(dim(generate_latent_values(n_1, p_1, p_1, s_1)),
                 c(n+nrow(p_1), 1))
    expect_equal(dim(generate_latent_values(n_2, p_2, p_2, s_2)),
                 c(n+nrow(p_2), 2))
})


sim_1 <- simulate_network(20, k=1)
sim_2 <- simulate_network(20, k=2)
sim_3 <- simulate_network(20, k=2, periods=2)

test_that("simulate_network produces the correct structure", {
    expect_equal(dim(positions(sim_1)), c(20, 1))
    expect_null(trajectories(sim_1))
    expect_equal(dim(positions(sim_2)), c(20, 2))
    expect_null(trajectories(sim_2))
    expect_equal(dim(trajectories(sim_3)), c(20, 2))
})

test_that("simulate_network produces warnings and errors", {
    ## n cannot be less than k+1
    expect_error(simulate_network(10, k=10))
    expect_error(simulate_network(10, k=2, ref_pos=matrix(0)))
    expect_error(simulate_network(10, k=2, mean_pos=matrix(0)))
    expect_error(simulate_network(10, k=2, groups=2,
                                  sigma_pos=list(matrix(0))))
    expect_warning(simulate_network(10, k=2,
                                    ref_pos=matrix(0, ncol=2, nrow=3)))
})
