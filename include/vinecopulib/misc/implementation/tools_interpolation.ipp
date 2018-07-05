// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <stdexcept>
#include <iostream>

namespace vinecopulib {

namespace tools_interpolation {
//! Constructor
//!
//! @param grid_points an ascending sequence of grid_points; used in both
//! dimensions.
//! @param values a dxd matrix of copula density values evaluated at
//! (grid_points_i, grid_points_j).
//! @param norm_times how many times the normalization routine should run.
inline InterpolationGrid::InterpolationGrid(const Eigen::VectorXd &grid_points,
                                            const Eigen::MatrixXd &values,
                                            int norm_times)
{
    if (values.cols() != values.rows()) {
        throw std::runtime_error("values must be a quadratic matrix");
    }
    if (grid_points.size() != values.rows()) {
        throw std::runtime_error(
            "number of grid_points must equal dimension of values");
    }

    grid_points_ = grid_points;
    values_ = values;
    normalize_margins(norm_times);
}

inline Eigen::MatrixXd InterpolationGrid::get_values() const
{
    return values_;
}

inline void InterpolationGrid::set_values(const Eigen::MatrixXd &values,
                                          int norm_times)
{
    if (values.size() != values_.size()) {
        if (values.rows() != values_.rows()) {
            std::stringstream message;
            message <<
                    "values have has wrong number of rows; " <<
                    "expected: " << values_.rows() << ", " <<
                    "actual: " << values.rows() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        if (values.cols() != values_.cols()) {
            std::stringstream message;
            message <<
                    "values have wrong number of columns; " <<
                    "expected: " << values_.cols() << ", " <<
                    "actual: " << values.cols() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }

    values_ = values;
    normalize_margins(norm_times);
}

inline void InterpolationGrid::flip()
{
    values_.transposeInPlace();
}

//! renormalizes the estimate to uniform margins
//!
//! @param times how many times the normalization routine should run.
inline void InterpolationGrid::normalize_margins(int times)
{
    size_t m = grid_points_.size();
    for (int k = 0; k < times; ++k) {
        for (size_t i = 0; i < m; ++i) {
            values_.row(i) /= int_on_grid(1.0, values_.row(i), grid_points_);
        }
        for (size_t j = 0; j < m; ++j) {
            values_.col(j) /= int_on_grid(1.0, values_.col(j), grid_points_);
        }
    }
}

inline Eigen::Matrix<ptrdiff_t, 1, 2> InterpolationGrid::get_indices(
    double x0, double x1)
{
    Eigen::Matrix<ptrdiff_t, 1, 2> out;
    out << 0, 0;
    bool found_i = false;
    bool found_j = false;
    for (ptrdiff_t k = 1; k < (grid_points_.size() - 1); ++k) {
        if ((x0 >= grid_points_(k))) {
            out(0) = k;
        } else {
            found_i = true;
        }
        if ((x1 >= grid_points_(k))) {
            out(1) = k;
        } else {
            found_j = true;
        }
        if (found_i & found_j) {
            break;
        }
    }
    return out;
}

inline double
bilinear_interpolation(double z11, double z12, double z21, double z22,
                       double x1, double x2, double y1, double y2,
                       double x, double y)
{
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return (z11 * x2x * y2y +
            z21 * xx1 * y2y +
            z12 * x2x * yy1 +
            z22 * xx1 * yy1) / (x2x1 * y2y1);
}

//! Interpolation in two dimensions
//!
//! @param x mx2 matrix of evaluation points.
inline Eigen::VectorXd
InterpolationGrid::interpolate(const Eigen::MatrixXd &x)
{

    auto f = [this](double x0, double x1) {

        auto indices = this->get_indices(x0, x1);
        return fmax(bilinear_interpolation(
            this->values_(indices(0), indices(1)),
            this->values_(indices(0), indices(1) + 1),
            this->values_(indices(0) + 1, indices(1)),
            this->values_(indices(0) + 1, indices(1) + 1),
            this->grid_points_(indices(0)),
            this->grid_points_(indices(0) + 1),
            this->grid_points_(indices(1)),
            this->grid_points_(indices(1) + 1),
            x0, x1), 1e-15);
    };

    return tools_eigen::binaryExpr_or_nan(x, f);
}

//! Integrate the grid along one axis
//!
//! @param u mx2 matrix of evaluation points
//! @param cond_var either 1 or 2; the axis considered fixed.
//!
inline Eigen::VectorXd
InterpolationGrid::integrate_1d(const Eigen::MatrixXd &u,
                                 size_t cond_var)
{
    ptrdiff_t m = grid_points_.size();
    Eigen::VectorXd tmpvals(m);
    Eigen::MatrixXd tmpgrid(m, 2);

    auto f = [this, m, cond_var, &tmpvals, &tmpgrid](double u1, double u2) {
        double upr = 0.0;
        double tmpint = 0.0, int1;
        if (cond_var == 1) {
            upr = u2;
            tmpgrid.col(0) = Eigen::VectorXd::Constant(m, u1);
            tmpgrid.col(1) = grid_points_;
        } else if (cond_var == 2) {
            upr = u1;
            tmpgrid.col(0) = grid_points_;
            tmpgrid.col(1) = Eigen::VectorXd::Constant(m, u2);
        }
        tmpvals = interpolate(tmpgrid).array().max(1e-4);
        tmpint = int_on_grid(upr, tmpvals, grid_points_);
        int1 = int_on_grid(1.0, tmpvals, grid_points_);

        return fmin(fmax(tmpint / int1, 1e-10), 1 - 1e-10);
    };

    return tools_eigen::binaryExpr_or_nan(u, f);
}

//! Integrate the grid along the two axis
//!
//! @param u mx2 matrix of evaluation points
//!
inline Eigen::VectorXd
InterpolationGrid::integrate_2d(const Eigen::MatrixXd &u)
{
    ptrdiff_t m = grid_points_.size();
    Eigen::VectorXd tmpvals(m), tmpvals2(m);
    Eigen::MatrixXd tmpgrid(m, 2);
    tmpgrid.col(1) = grid_points_;

    auto f = [this, m, &tmpvals, &tmpvals2, &tmpgrid](double u1, double u2) {
        double upr, tmpint, tmpint1;
        upr = u2;
        for (ptrdiff_t k = 0; k < m; ++k) {
            tmpgrid.col(0) = Eigen::VectorXd::Constant(m, grid_points_(k));
            tmpvals = interpolate(tmpgrid);
            tmpint = int_on_grid(upr, tmpvals, grid_points_);
            tmpvals2(k) = tmpint;
        }
        upr = u1;
        tmpint = int_on_grid(upr, tmpvals2, grid_points_);
        tmpint1 = int_on_grid(1.0, tmpvals2, grid_points_);
        return fmin(fmax(tmpint / tmpint1, 1e-10), 1 - 1e-10);
    };

    return tools_eigen::binaryExpr_or_nan(u, f);
}

// ---------------- Utility functions for spline interpolation ----------------

//! Evaluate a cubic polynomial
//!
//! @param x evaluation point.
//! @param a polynomial coefficients
inline double
InterpolationGrid::cubic_poly(const double &x, const Eigen::VectorXd &a)
{
    double x2 = x * x;
    double x3 = x2 * x;
    return a(0) + a(1) * x + a(2) * x2 + a(3) * x3;
}

//! Indefinite integral of a cubic polynomial
//!
//! @param x evaluation point.
//! @param a polynomial coefficients.
inline double InterpolationGrid::cubic_indef_integral(const double &x,
                                                      const Eigen::VectorXd &a)
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    return a(0) * x + a(1) / 2.0 * x2 + a(2) / 3.0 * x3 + a(3) / 4.0 * x4;
}

//! Definite integral of a cubic polynomial
//!
//! @param lower lower limit of the integral.
//! @param upper upper limit of the integral.
//! @param a polynomial coefficients.
inline double InterpolationGrid::cubic_integral(const double &lower,
                                                const double &upper,
                                                const Eigen::VectorXd &a)
{
    return cubic_indef_integral(upper, a) - cubic_indef_integral(lower, a);
}

//! Calculate coefficients for cubic intrpolation spline
//!
//! @param vals length 4 vector of function values.
//! @param grid length 4 vector of grid points.
inline Eigen::VectorXd
InterpolationGrid::find_coefs(const Eigen::VectorXd &vals,
                              const Eigen::VectorXd &grid)
{
    Eigen::VectorXd a(4);

    double dt0 = grid(1) - grid(0);
    double dt1 = grid(2) - grid(1);
    double dt2 = grid(3) - grid(2);

    /* check for repeated points (important for boundaries) */
    if (dt1 < 1e-4)
        dt1 = 1.0;
    if (dt0 < 1e-4)
        dt0 = dt1;
    if (dt2 < 1e-4)
        dt2 = dt1;

    // compute tangents when parameterized in (t1,t2)
    double dx1 = (vals(1) - vals(0)) / dt0;
    dx1 -= (vals(2) - vals(0)) / (dt0 + dt1);
    dx1 += (vals(2) - vals(1)) / dt1;
    double dx2 = (vals(2) - vals(1)) / dt1;
    dx2 -= (vals(3) - vals(1)) / (dt1 + dt2);
    dx2 += (vals(3) - vals(2)) / dt2;

    // rescale tangents for parametrization in (0,1)
    dx1 *= dt1;
    dx2 *= dt1;

    // compute coefficents
    a(0) = vals(1);
    a(1) = dx1;
    a(2) = -3 * vals(1) + 3 * vals(2) - 2 * dx1 - dx2;
    a(3) = 2 * vals(1) - 2 * vals(2) + dx1 + dx2;

    return a;
}

//! Interpolate on 4 points
//!
//! @param x evaluation point.
//! @param vals length 4 vector of function values.
//! @param grid length 4 vector of grid points.
inline double InterpolationGrid::interp_on_grid(const double &x,
                                                const Eigen::VectorXd &vals,
                                                const Eigen::VectorXd &grid)
{
    Eigen::VectorXd a = find_coefs(vals, grid);
    double xev = fmax((x - grid(1)), 0) / (grid(2) - grid(1));
    return cubic_poly(xev, a);
}


// ---------------- Utility functions for integration ----------------


//! Integrate a spline interpolant
//!
//! @param upr upper limit of integration (lower is 0).
//! @param vals vector of values to be interpolated and integrated.
//! @param grid vector of grid points on which vals has been computed.
//!
//! @return Integral of interpolation spline defined by (vals, grid).
inline double InterpolationGrid::int_on_grid(const double &upr,
                                             const Eigen::VectorXd &vals,
                                             const Eigen::VectorXd &grid)
{
    ptrdiff_t m = grid.size();
    Eigen::VectorXd tmpvals(4), tmpgrid(4), tmpa(4), a(4);
    double uprnew, newint;

    double tmpint = 0.0;

    if (upr > grid(0)) {
        // go up the grid and integrate
        for (ptrdiff_t k = 0; k < m - 1; ++k) {
            // stop loop if fully integrated
            if (upr < grid(k))
                break;

            // select length 4 subvectors and calculate spline coefficients
            tmpvals(0) = vals(std::max(k - 1, static_cast<ptrdiff_t>(0)));
            tmpvals(1) = vals(k);
            tmpvals(2) = vals(k + 1);
            tmpvals(3) = vals(std::min(k + 2, m - 1));

            tmpgrid(0) = grid(std::max(k - 1, static_cast<ptrdiff_t>(0)));
            tmpgrid(1) = grid(k);
            tmpgrid(2) = grid(k + 1);
            tmpgrid(3) = grid(std::min(k + 2, m - 1));

            tmpa = find_coefs(tmpvals, tmpgrid);

            // don't integrate over full cell if upr is in interior
            uprnew = (upr - grid(k)) / (grid(k + 1) - grid(k));
            newint = cubic_integral(0.0, fmin(1.0, uprnew), tmpa);
            tmpint += newint * (grid(k + 1) - grid(k));
        }
    }

    return tmpint;
}
}

}

//inline Eigen::VectorXd
//bilinear_interpolation(const Eigen::VectorXd& q11,
//                       const Eigen::VectorXd& q12,
//                       const Eigen::VectorXd& q21,
//                       const Eigen::VectorXd& q22,
//                       const Eigen::VectorXd& x1,
//                       const Eigen::VectorXd& x2,
//                       const Eigen::VectorXd& y1,
//                       const Eigen::VectorXd& y2,
//                       const Eigen::VectorXd& x,
//                       const Eigen::VectorXd& y)
//{
//    size_t n = q11.size();
//    Eigen::VectorXd x2x1(n), y2y1(n), x2x(n), y2y, yy1(n), xx1(n);
//    x2x1 = x2 - x1;
//    y2y1 = y2 - y1;
//    x2x = x2 - x;
//    y2y = y2 - y;
//    yy1 = y - y1;
//    xx1 = x - x1;
//    return (q11.cwiseProduct(x2x).cwiseProduct(y2y) +
//        q21.cwiseProduct(xx1).cwiseProduct(y2y) +
//        q12.cwiseProduct(x2x).cwiseProduct(yy1) +
//        q22.cwiseProduct(xx1).cwiseProduct(yy1)).cwiseQuotient(x2x1.cwiseProduct(y2y1));
//}
//
//inline double
//bilinear_integration(double z11, double z12, double z21, double z22,
//                     double x1, double x2, double y1, double y2)
//{
//    return (z11 + z21 + z12 + z22) * (x2 - x1) * (y2 - y1) / 4;
//}
//
//inline Eigen::VectorXd
//InterpolationGrid::integrate_2d(const Eigen::MatrixXd &u) {
//
//    ptrdiff_t m = grid_points_.size();
//
//    auto f = [this, m](double u1, double u2) {
//
//        double res = 0;
//        double x1, x2, y1, y2, z11, z12, z21, z22;
//        for (ptrdiff_t k1 = 0; k1 < (m - 1); ++k1) {
//
//            x1 = this->grid_points_(k1);
//            x2 = this->grid_points_(k1 + 1);
//
//            for (ptrdiff_t k2 = 0; k2 < (m - 1); ++k2) {
//
//                y1 = this->grid_points_(k2);
//                y2 = this->grid_points_(k2 + 1);
//                z11 = this->values_(k1, k2);
//
//                if ((u1 >= x2) || (u2 >= y2)) {
//                    z12 = bilinear_interpolation(this->values_(k1, k2),
//                                                 this->values_(k1, k2 + 1),
//                                                 this->values_(k1 + 1, k2 ),
//                                                 this->values_(k1 + 1, k2 + 1),
//                                                 x1, x2, y1, y2, x1, u2);
//                    z21 = bilinear_interpolation(this->values_(k1, k2),
//                                                 this->values_(k1, k2 + 1),
//                                                 this->values_(k1 + 1, k2 ),
//                                                 this->values_(k1 + 1, k2 + 1),
//                                                 x1, x2, y1, y2, u1, y1);
//                    z22 = bilinear_interpolation(this->values_(k1, k2),
//                                                 this->values_(k1, k2 + 1),
//                                                 this->values_(k1 + 1, k2 ),
//                                                 this->values_(k1 + 1, k2 + 1),
//                                                 x1, x2, y1, y2, u1, y2);
//                } else {
//                    z12 = this->values_(k1, k2 + 1);
//                    z21 = this->values_(k1 + 1, k2 );
//                    z22 = this->values_(k1 + 1, k2 + 1);
//                }
//
//                res += bilinear_integration(z11, z12, z21, z22, x1, x2, y1, y2);
//                if (u2 >= x2) {
//                    break;
//                }
//            }
//            if ((u1 >= x1)) {
//                break;
//            }
//        }
//
//        return fmin(fmax(res, 1e-10), 1 - 1e-10);
//    };
//
//    return tools_eigen::binaryExpr_or_nan(u, f);
//}

//inline Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic> InterpolationGrid::get_indices2(const Eigen::MatrixXd &x)
//{
//    Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic> indices(x.rows(), x.cols());
//    indices.fill(0);
//
//    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> found(x.rows(), x.cols());
//    found.fill(false);
//
//    for (ptrdiff_t k = 1; k < (grid_points_.size() - 1); ++k) {
//
//        found = (x.array() < grid_points_(k)).select(true, found);
//        indices = (found.array() != true).select(k, indices);
//
//        if (found.all())
//            break;
//    }
//
//    return indices;
//}
//
//inline Eigen::Matrix<ptrdiff_t, 1, 2> InterpolationGrid::get_indices(
//    double x0, double x1)
//{
//
//    Eigen::Matrix<ptrdiff_t, 1, 2> out;
//    out << 0, 0;
//    out(0) = tools_eigen::interpolation_search(grid_points_, x0);
//    out(1) = tools_eigen::interpolation_search(grid_points_, x1);
//
//    return out;
//}