## helper function
get_predvars <- function(mod, key_attr) {
    out <- as.list(attr(terms(mod$model), "predvars")[[3]])[key_attr]
    if (any(sapply(out, is.null)))
        stop("Found not matched key attribute.")
    out
}
get_attr <- function(x, key_attr) {
    out <- attributes(x)[key_attr]
    if (any(sapply(out, is.null)))
        stop("Found not matched key attribute.")
    out
}

## simulated data
n <- 1e2
x <- rnorm(n)
y <- x + rnorm(n, sd = 0.1)
new_x <- runif(2 * n, min(x), max(x))

## bSpline()
mod <- lm(y ~ bsp(x, df = 6))
key_attr <- c("degree", "knots", "Boundary.knots", "intercept",
              "periodic", "derivs", "integral")
expect_equal(
    get_attr(bSpline(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

new_mat <- predict(bSpline(x, df = 6), new_x)
pred1 <- predict(mod, data.frame(x = new_x))
pred2 <- coef(mod)[1L] + as.numeric(new_mat %*% coef(mod)[- 1L])
expect_equivalent(pred1, pred2)

## design matrix
X <- bsp(x, df = 6)
mod <- lm(y ~ X)
expect_error(get_predvars(mod, key_attr), "not matched key attribute")

## naturalSpline()
mod <- lm(y ~ naturalSpline(x, df = 6))
key_attr <- c("knots", "Boundary.knots", "intercept")
expect_equal(
    get_attr(nsp(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

## design matrix
X <- nsp(x, df = 6)
mod <- lm(y ~ X)
expect_error(get_predvars(mod, key_attr), "not matched key attribute")

## mSpline()
mod <- lm(y ~ mSpline(x, df = 6))
key_attr <- c("degree", "knots", "Boundary.knots", "intercept",
              "periodic", "derivs", "integral")
expect_equal(
    get_attr(mSpline(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

## iSpline()
mod <- lm(y ~ iSpline(x, df = 6))
key_attr <- c("degree", "knots", "Boundary.knots", "intercept", "derivs")
expect_equal(
    get_attr(iSpline(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

## cSpline()
mod <- lm(y ~ cSpline(x, df = 6))
key_attr <- c("degree", "knots", "Boundary.knots", "intercept",
              "scale", "derivs")
expect_equal(
    get_attr(cSpline(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

## dbs()
mod <- lm(y ~ dbs(x, df = 6))
key_attr <- c("degree", "knots", "Boundary.knots", "intercept", "derivs")
expect_equal(
    get_attr(dbs(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

## ibs()
mod <- lm(y ~ ibs(x, df = 6))
key_attr <- c("degree", "knots", "Boundary.knots", "intercept")
expect_equal(
    get_attr(ibs(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)

## bernsteinPoly()
mod <- lm(y ~ bernsteinPoly(x, df = 6))
key_attr <- c("degree", "Boundary.knots", "intercept", "derivs", "integral")
expect_equal(
    get_attr(bernsteinPoly(x, df = 6), key_attr),
    get_predvars(mod, key_attr)
)
