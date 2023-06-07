library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## naturalSpline()
nsMat0 <- naturalSpline(x, knots = knots, intercept = TRUE)
nsMat1 <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
nsMat2 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
nsMat3 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 2)

op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
plot(nsMat0, ylab = "basis")
plot(nsMat1, ylab = "integral")
plot(nsMat2, ylab = "1st derivative")
plot(nsMat3, ylab = "2nd derivative")
par(op) # reset to previous plotting settings

## nsk()
nskMat <- nsk(x, knots = knots, intercept = TRUE)
plot(nskMat, ylab = "nsk()", mark_knots = "all")
abline(h = 1, col = "red", lty = 3)

## use the deriv method
all.equal(nsMat0, deriv(nsMat1), check.attributes = FALSE)
all.equal(nsMat2, deriv(nsMat0))
all.equal(nsMat3, deriv(nsMat2))
all.equal(nsMat3, deriv(nsMat0, 2))

## a linear model example
fit1 <- lm(weight ~ -1 + nsk(height, df = 4, intercept = TRUE), data = women)
fit2 <- lm(weight ~ nsk(height, df = 3), data = women)

## the knots (same for both fits)
knots <- unlist(attributes(fit1$model[[2]])[c('Boundary.knots', 'knots')])

## predictions at the knot points
predict(fit1, data.frame(height = sort(unname(knots))))
unname(coef(fit1)) # equal to the coefficient estimates

## different interpretation when "intercept = FALSE"
unname(coef(fit1)[-1] - coef(fit1)[1]) # differences: yhat[2:4] - yhat[1]
unname(coef(fit2))[-1]                 # ditto
