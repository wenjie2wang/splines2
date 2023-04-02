## check essentials
expect_eqt <- function(current, target, ...)
{
    tinytest::expect_equivalent(current = unclass(current),
                                target = unclass(target), ...)
}
