//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2025
//
// This file is part of the R package splines2.
//
// The R package splines2 is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package splines2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#ifndef SPLINES2_UTILS_H
#define SPLINES2_UTILS_H

#include <algorithm>            // std::max, std::set_union, etc.
#include <cmath>                // std::pow and std::sqrt, etc.
#include <limits>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

// #include <armadillo>
#include <RcppArmadillo.h>

namespace splines2 {

    // compare double-precision numbers for almost equality
    inline bool isAlmostEqual(double A, double B = 0.0)
    {
        const double MaxRelDiff { std::numeric_limits<double>::epsilon() };
        // compute the difference.
        double diff = std::abs(A - B);
        A = std::abs(A);
        B = std::abs(B);
        // Find the largest
        double largest = (B > A) ? B : A;
        if (diff <= largest * MaxRelDiff) {
            return true;
        }
        return false;
    }

    template <typename T>
    inline bool is_approx_equal(const T& A, const T& B) {
        const double MaxRelDiff { std::numeric_limits<double>::epsilon() };
        return arma::approx_equal(A, B, "reldiff", MaxRelDiff);
    }

    template <typename T>
    // function checking if there exists any duplicates
    inline bool any_duplicated(const T& x)
    {
        std::unordered_set<double> seen;
        bool res {false};
        for (size_t i {0}; i < x.n_rows; ++i) {
            res = ! seen.insert(x(i)).second;
            if (res) break;
        }
        return res;
    }

    // cumulative sum in possibly reverse order>
    template <typename T>
    inline T rev_cum_sum(const T& x)
    {
        const unsigned long int n_x {x.n_elem};
        T res { arma::zeros<T>(n_x) };
        double tmp {0.0};
        for (size_t i {1}; i <= n_x; ++i) {
            tmp += x[n_x - i];
            res[n_x - i] = tmp;
        }
        return res;
    }

    template <typename T>
    inline T get_inside_x(const T& x, const arma::vec& boundary_knots)
    {
        std::vector<double> res;
        for (size_t i {0}; i < x.n_elem; ++i) {
            if (x(i) >= boundary_knots(0) && x(i) <= boundary_knots(1)) {
                res.push_back(x(i));
            }
        }
        return T(res);
    }

    // quantile function
    // reference: Hyndman and Fan (1996)
    inline arma::vec quantile(const arma::vec& x,
                              const arma::vec& probs,
                              const unsigned int type = 7)
    {
        double alpha { 0 }, beta { 0 };
        switch(type) {
            case 4:
                beta = 1.0;
                break;
            case 5:
                alpha = beta = 0.5;
                break;
            case 6:
                break;
            default:
            case 7:
                alpha = beta = 1.0;
                break;
            case 8:
                alpha = beta = 1.0 / 3;
                break;
            case 9:
                alpha = beta = 3.0 / 8;
                break;
        }
        const long n { static_cast<long>(x.n_elem) };
        arma::vec inc_x { arma::sort(x) };
        arma::vec res { arma::zeros(probs.n_elem) };
        double fuzz { std::numeric_limits<double>::epsilon() };
        for (size_t i {0}; i < probs.n_elem; ++i) {
            // n * p + m
            double nppm { alpha + probs(i) * (n + 1 - alpha - beta) };
            double j { std::floor(nppm + fuzz) };
            double h { nppm - j };
            long lj { static_cast<long>(j) };
            if (lj == 0) {
                res(i) = x.min();
                continue;
            }
            if (lj >= n) {
                res(i) = x.max();
                continue;
            }
            res(i) = (1 - h) * inc_x(lj - 1) + h * inc_x(lj);
        }
        return res;
    }

    // x of length inside (start, end) from seq(start, end, length.out + 2)
    inline arma::vec linspace_inside(const double start,
                                     const double end,
                                     const unsigned int length)
    {
        arma::vec out { arma::linspace(start, end, length + 2) };
        return out.subvec(1, length);
    }

    // inline handy functions
    inline arma::vec mat2vec(const arma::mat& x) {
        return arma::conv_to<arma::vec>::from(x);
    }
    inline arma::rowvec mat2rowvec(const arma::mat& x) {
        return arma::conv_to<arma::rowvec>::from(x);
    }
    inline arma::vec num2vec(const double x) {
        arma::vec out { arma::zeros(1) };
        out(0) = x;
        return out;
    }

    // convert arma vector type to Rcpp vector type
    template <typename T>
    inline Rcpp::NumericVector arma2rvec(const T& x) {
        return Rcpp::NumericVector(x.begin(), x.end());
    }
    template <typename T>
    inline Rcpp::IntegerVector arma2ivec(const T& x) {
        return Rcpp::IntegerVector(x.begin(), x.end());
    }
    // convert Rcpp::NumericVector to arma::colvec
    template <typename T>
    inline arma::vec rvec2arma(const T& x) {
        return arma::vec(x.begin(), x.size(), false);
    }
    // convert arma matrix type to Rcpp matrix type
    inline Rcpp::NumericMatrix arma2rmat(const arma::mat& x) {
        return Rcpp::NumericMatrix(x.n_rows, x.n_cols, x.begin());
    }
    // convert Rcpp matrix to arma matrix
    // cannot add const to x
    inline arma::mat rmat2arma(Rcpp::NumericMatrix& x) {
        return arma::mat(x.begin(), x.nrow(), x.ncol(), false);
    }

    // function to remove the first column of a matrix
    template <typename T>
    inline T mat_wo_col1(const T& x)
    {
        const arma::uword x_ncol { x.n_cols };
        if (x_ncol > 1) {
            return x.tail_cols(x_ncol - 1);
        }
        // else
        throw std::range_error("No column left in the matrix.");
        return T();
    }

    // function to add zero columns to the end of a matrix
    inline arma::mat add_zero_cols(const arma::mat& x,
                                   const unsigned int n_cols = 1)
    {
        // create zero matrix
        arma::mat mat2 { arma::zeros(x.n_rows, n_cols) };
        return arma::join_rows(x, mat2);
    }

    // create a character vector ("1", ...) of a given length
    inline Rcpp::CharacterVector char_seq_len(const unsigned int n) {
        Rcpp::CharacterVector out { n };
        for (size_t i {0}; i < n; ++i) {
            out[i] = std::to_string(i + 1);
        }
        return out;
    }

}

#endif
