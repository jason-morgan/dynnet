dapply <- function(n, fn, ...)
{
    output <- vector("list", length=n*(n-1)/2)
    idx <- 1

    for (.J in 1:(n-1)) {
        for (.I in (.J+1):n) {
            environment(fn) <- environment()
            output[[idx]] <- fn(...)
            idx <- idx + 1
        }
    }

    output
}

fordyads <- function(n, fn, ...)
{
    for (.J in 1:(n-1)) {
        for (.I in (.J+1):n) {
            environment(fn) <- environment()
            fn(...)
        }
    }

    invisible()
}
