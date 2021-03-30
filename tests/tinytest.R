if (requireNamespace("tinytest", quietly = TRUE) &&
    utils::packageVersion("tinytest") >= "1.2.2") {

    ## Set a seed to make the test deterministic
    set.seed(808)
    is_at_home <- function() {
        identical(tolower(Sys.getenv("TT_AT_HOME")),"true")
    }

    ## only run the following tests if at home (not at cran)
    tinytest::test_package("splines2", ncpu = NULL,
                           test_dir = "rcpp-tests",
                           at_home = is_at_home())
    ## run anyway
    tinytest::test_package("splines2", ncpu = NULL,
                           side_effects = TRUE)
}
