if (requireNamespace("tinytest", quietly = TRUE) &&
    utils::packageVersion("tinytest") >= "1.2.2") {
    set.seed(808)
    is_at_home <- function() {
        identical(tolower(Sys.getenv("TT_AT_HOME")), "true")
    }
    ## only run the following tests if "at home" (not at cran)
    tinytest::test_package("splines2", ncpu = NULL,
                           testdir = "rcpp-tests",
                           at_home = is_at_home())
    ## always run the following checks
    tinytest::test_package("splines2", ncpu = NULL,
                           side_effects = TRUE)
}
