context("Network Constructors and Accessors")

x1 <- matrix(c(0, 1, 0,
               1, 0, 1,
               0, 1, 0), ncol=3, byrow=TRUE)

x2 <- matrix(c(0.0, 1.1, 1.2,
               1.1, 0.0, 3.0,
               0.5, 0.3, 0.0), ncol=3, byrow=TRUE)

net_b  <- dynnet_adjacency(x1)
net_b2 <- dynnet_adjacency(list(x, x), mode="directed", weighted=NULL)

net_w  <- dynnet_adjacency(x, mode="directed", weighted=TRUE)
net_w2 <- dynnet_adjacency(list(x, x), mode="directed", weighted=TRUE)

test_that("dynnet object has correct structure", {
    expect_equal(net_b$nodes, 3)
})

test_that("get_adjacency returns the adjacency matrices", {
    expect_equal(get_adjacency(net_w, period=1)[2,1], x2[2,1])
})
