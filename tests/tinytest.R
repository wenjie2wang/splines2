if (requireNamespace("tinytest", quietly = TRUE) &&
    utils::packageVersion("tinytest") >= "1.0.0") {

    ## Set a seed to make the test deterministic
    set.seed(808)

    tinytest::test_package("splines2", ncpu = getOption("Ncpus", 1),
                           side_effects = TRUE)
}
