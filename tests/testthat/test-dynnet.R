context("dynnet object constructors and accessors")

x <- matrix(c(0.0, 1.1, 1.2,
              1.1, 0.0, 3.0,
              0.5, 0.3, 0.0), ncol=3, byrow=TRUE)

net_b  <- dynnet_adjacency(x, directed=TRUE, weighted=FALSE)
net_b2 <- dynnet_adjacency(list(x, x), directed=TRUE, weighted=FALSE)
net_w  <- dynnet_adjacency(x, directed=TRUE, weighted=TRUE)
net_w2 <- dynnet_adjacency(list(x, x), directed=TRUE, weighted=TRUE)

test_that("get_adjacency returns the adjacency matrices", {
    expect_equal(get_adjacency(net_w)[[1]], x)
    expect_equal(get_adjacency(net_w, 1), x)
})

test_that("dynnet_adjacency is producing binary and weighted networks", {
    expect_equal(get_adjacency(net_b, 1)[1,1],  0)    
    expect_equal(get_adjacency(net_b, 1)[2,1],  1)
    expect_equal(get_adjacency(net_w, 1)[2,1],  1.1)    
    expect_equal(get_adjacency(net_b2, 2)[2,1], 1)
    expect_equal(get_adjacency(net_w2), list(x, x))
    expect_equal(get_adjacency(net_w2, 1:2), list(x, x))        
})
