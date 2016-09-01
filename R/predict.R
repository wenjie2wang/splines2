##'
predict.ibs <- function(object, newx, ...) {
    if (missing(newx)) return(object)
    a <- c(list(x = newx),
          attributes(object)[c("degree", "knots", "Boundary.knots",
                               "intercept")])
    do.call("ibs", a)
}
