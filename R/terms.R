absdiff <- function(graph, vattr)
{
    x <- vertex_attr(graph, vattr)
    if (is.factor(x))
        stop("variable ", vattr,
             " is a factor. absdiff supports real values only")

    a <- outer(x, x, function(x1, x2) abs(x1 - x2))
    a[lower.tri(a)]
}
