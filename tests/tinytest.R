if (requireNamespace("tinytest", quietly = TRUE) &&
    utils::packageVersion("tinytest") >= "1.2.2") {

    ## Set a seed to make the test deterministic
    set.seed(808)

    ## only run the following tests if not at cron
    tinytest::test_package("splines2", ncpu = NULL,
                           test_dir = "rcpp-tests",
                           at_home = tinytest::at_home())
    ## run anyway
    tinytest::test_package("splines2", ncpu = NULL,
                           side_effects = TRUE)
}
