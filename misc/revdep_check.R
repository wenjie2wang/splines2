revdep_dir <- "revdep"
res <- tools::check_packages_in_dir(revdep_dir,
                                    clean = FALSE,
                                    reverse = list(recursively = FALSE))
summary(res)
