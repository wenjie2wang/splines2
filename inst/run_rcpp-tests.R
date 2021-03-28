#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (! is.na(args)) setwd(args[1L])

tinytest::run_test_dir(dir = "rcpp-tests")
