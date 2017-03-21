// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "interpolation_grid.hpp"
#include "tools_stats.hpp"
#include <exception>

namespace vinecopulib
{
    //! Constructor
    //!
    //! @param grid_points an ascending sequence of grid_points; used in both
    //! dimensions.
    //! @param values a dxd matrix of copula density values evaluated at
    //! (grid_points_i, grid_points_j).
    InterpolationGrid::InterpolationGrid(
        const Eigen::VectorXd& grid_points, 
        const Eigen::MatrixXd& values
    )
    {
        if (values.cols() != values.rows()) {
            throw std::runtime_error(
                    "values must be a quadratic matrix"
            );
        }
        if (grid_points.size() != values.rows()) {
            throw std::runtime_error(
                    "number of grid_points must equal dimension of values"
            );
        }

        grid_points_ = grid_points;
        values_ = values;
    }

    void InterpolationGrid::flip()
    {
        values_.transposeInPlace();
    }

    //! Interpolation in two dimensions
    //!
    //! @param x mx2 matrix of evaluation points.
    Eigen::VectorXd InterpolationGrid::interpolate(const Eigen::MatrixXd& x)
    {
        int N = x.rows();
        int m = grid_points_.size();
        Eigen::VectorXd y(4), out(N), a(4), tmpgrid(4), tmpvals(4);
        int i = 0;
        int j = 0;
        int i0, i3;

        for (int n = 0; n < N; ++n) {
            // find cell
            bool found_i = false;
            bool found_j = false;
            for (int k = 1; k < (m-1); ++k) {
                if ((x(n, 0) >= grid_points_(k))) { 
                    i = k;
                } else {
                    found_i = true;
                }
                if ((x(n, 1) >= grid_points_(k))) {
                    j = k;
                } else {
                    found_j = true;
                }
                if (found_i & found_j) {
                    break;
                }
            }

            // construct grid for first direction
            i0 = std::max(i-1, 0);
            i3 = std::min(i+2, m-1);
            tmpgrid(0) = grid_points_(i0);
            tmpgrid(1) = grid_points_(i);
            tmpgrid(2) = grid_points_(i+1);
            tmpgrid(3) = grid_points_(i3);

            // interpolate in one direction (four times)
            for (int s = 0; s < 4; ++s) {
                i0 = std::max(i-1, 0);
                i3 = std::min(i+2, m-1);
                int jj = std::min(m-1, j-1+s);
                jj = std::max(0, jj);

                tmpvals(0) = values_(i0,  jj);
                tmpvals(1) = values_(i,   jj);
                tmpvals(2) = values_(i+1, jj);
                tmpvals(3) = values_(i3,  jj);

                y(s) = interp_on_grid(x(n, 0), tmpvals, tmpgrid);
                y(s) = fmax(y(s), 0.0);
            }

            // use these four points to interpolate in the remaining direction#
            i0 = std::max(j-1, 0);
            i3 = std::min(j+2, m-1);
            tmpgrid(0) = grid_points_(i0);
            tmpgrid(1) = grid_points_(j);
            tmpgrid(2) = grid_points_(j+1);
            tmpgrid(3) = grid_points_(i3);

            out(n) = interp_on_grid(x(n, 1), y, tmpgrid);
            out(n) = fmax(out(n), 1e-15);
        }

        return out;
    }

    //! Integrate the grid along one axis
    //!
    //! @param u mx2 matrix of evaluation points
    //! @param cond_var either 1 or 2; the axis considered fixed.
    //!
    Eigen::VectorXd InterpolationGrid::intergrate_1d(
        const Eigen::MatrixXd& u, 
        const int& cond_var
    )
    {
        int n = u.rows();
        int m = grid_points_.size();
        Eigen::VectorXd tmpvals(m), out(n), tmpa(4), tmpb(4);
        Eigen::MatrixXd tmpgrid(m, 2);
        double upr = 0.0;
        double tmpint, int1;
        tmpint = 0.0;

        for (int i = 0; i < n; ++i) {
            if (cond_var == 1) {
                upr = u(i, 1);
                tmpgrid.col(0) = Eigen::VectorXd::Constant(m, u(i, 0));
                tmpgrid.col(1) = grid_points_;
            } else if (cond_var == 2) {
                upr = u(i, 0);
                tmpgrid.col(0) = grid_points_;
                tmpgrid.col(1) = Eigen::VectorXd::Constant(m, u(i, 1));
            }
            tmpvals = interpolate(tmpgrid);
            tmpint = int_on_grid(upr, tmpvals, grid_points_);
            int1 =  int_on_grid(1.0, tmpvals, grid_points_);
            out(i) = tmpint/int1;
            out(i) = fmax(out(i), 1e-10);
            out(i) = fmin(out(i), 1-1e-10);
        }

        return out;
    }

    //! Inverse of integral along one axis of the grid
    //!
    //! @param u mx2 matrix of evaluation points
    //! @param cond_var either 1 or 2; the axis considered fixed.
    //!
    Eigen::VectorXd InterpolationGrid::inv_intergrate_1d(
        const Eigen::MatrixXd& u,
        const int& cond_var
    )
    {
        Eigen::VectorXd out(u.rows());

        for (int i = 0; i < u.rows(); ++i) {
            int br = 0;
            double x0 = 0.0;
            double x1 = 1.0;

            // evaluation points at boundary
            Eigen::MatrixXd tmpu0(1, 2), tmpu1(1, 2);
            double q;
            if (cond_var == 1) {
                q = u(i, 1);
                tmpu0(0, 0) = u(i, 0);
                tmpu0(0, 1) = x0;
                tmpu1(0, 0) = u(i, 0);
                tmpu1(0, 1) = x1;

            } else if (cond_var == 2) {
                q = u(i, 0);
                tmpu0(0, 0) = x0;
                tmpu0(0, 1) = u(i, 1);
                tmpu1(0, 0) = x1;
                tmpu1(0, 1) = u(i, 1);
            } else {
                throw std::runtime_error(
                        "cond_var must be 1 or 2"
                );
            }

            // evaluate h-function at boundary points
            double ql = intergrate_1d(tmpu0, cond_var)(0);
            double qh = intergrate_1d(tmpu1, cond_var)(0);
            ql = ql - q;
            qh = qh - q;

            // check if already close enough (unless at boundary)
            double tol = 1e-10;
            double ans = 0.0;
            double val = 0.0;
            if ((::fabs(ql) < tol) && (q > 1e-9)) {
                ans = x0;
                br = 1;
            } else if ((::fabs(qh) < tol) && (q < 1-1e-9)) {
                ans = x1;
                br = 1;
            }

            // find inverse by bisection
            int maxit = 15;
            for (int it = 0; it < maxit; ++it) {
                // set new evaluation point
                ans = (x0 + x1) / 2.0;
                if (cond_var == 1) {
                    q = u(i, 1);
                    tmpu0(0, 0) = u(i, 0);
                    tmpu0(0, 1) = ans;

                } else if (cond_var == 2) {
                    q = u(i, 0);
                    tmpu0(0, 0) = ans;
                    tmpu0(0, 1) = u(i, 1);
                }

                double gap = intergrate_1d(tmpu0, cond_var)(0) - q;

                // find section for next iteration
                if (::fabs(gap) < 1e-9) {
                    if (q <= 9e-9) {
                        // go to upper section if q == 1e-10
                        x0 = ans;
                        ql = val;
                    } else if (q >= 1 - 9e-9)  {
                        // go to lower section if q == 1 - 1e-10
                        x1 = ans;
                        qh = val;
                    } else {
                        br = 1;
                    }
                } else if (val > 0.0) {
                    x1 = ans;
                    qh = gap;
                } else if (val < 0.0) {
                    x0 = ans;
                    ql = gap;
                }

                // stop if values are close enough
                if (::fabs(x0 - x1) <= tol) {
                    br = 1;
                }

                if (br == 1) {
                    break;
                }
            }

            out(i) = ans;
        }

        return out;
    }


// ---------------- Utility functions for spline interpolation ----------------

    //! Evaluate a cubic polynomial
    //!
    //! @param x evaluation point.
    //! @param a polynomial coefficients
    double InterpolationGrid::cubic_poly(
        const double& x, 
        const Eigen::VectorXd& a
    )
    {
        double x2 = x*x;
        double x3 = x2*x;
        return a(0) + a(1)*x + a(2)*x2 + a(3)*x3;
    }

    //! Indefinite integral of a cubic polynomial
    //!
    //! @param x evaluation point.
    //! @param a polynomial coefficients.
    double InterpolationGrid::cubic_indef_integral(
        const double& x,
        const Eigen::VectorXd& a
    )
    {
        double x2 = x*x;
        double x3 = x2*x;
        double x4 = x3*x;
        return a(0)*x + a(1)/2.0*x2 + a(2)/3.0*x3 + a(3)/4.0*x4;
    }

    //! Definite integral of a cubic polynomial
    //!
    //! @param lower lower limit of the integral.
    //! @param upper upper limit of the integral.
    //! @param a polynomial coefficients.
    double InterpolationGrid::cubic_integral(
        const double& lower, 
        const double& upper, 
        const Eigen::VectorXd& a
    )
    {
        return cubic_indef_integral(upper, a) - cubic_indef_integral(lower, a);
    }

    //! Numerically invert a cubic integral (with 0 as lower bound)
    //!
    //! @param q evaluation point (a 'quantile').
    //! @param a vector of polynomial coefficients.
    //!
    //! The inverse is found by the bisection method with a maximum of 20
    //! iterations.
    double InterpolationGrid::inv_cubic_integral(
        const double& q, 
        const Eigen::VectorXd& a
    )
    {
        double x0, x1, ql, qh, ans, val;
        ans = 0.0, val = 0.0; x0 = 0.0; x1 = 1.0;
        ql = 0.0;
        qh = cubic_integral(0.0, x1, a);
        int br = 0;
        double tol = ::fmax(1e-10 * (x1 - x0), 1e-10);

        // check if already close enough (or 1.0 is exceeded)
        ql = ql - q;
        qh = qh - q;
        if (::fabs(ql) <= tol) {
            ans = x0;
            br = 1;
        }
        if ((::fabs(qh) <= tol) | (qh < 0)) {
            ans = x1;
            br = 1;
        }

        for (int it = 0; it < 20; ++it) {
            ans = (x0 + x1) / 2.0;
            val = cubic_integral(0.0, ans, a);
            val = val - q;
            // stop if values become too close (avoid infinite loop)
            if (::fabs(val) <= tol)
                br = 1;
            if (::fabs(x0 - x1) <= tol)
                br = 1;
            // check which of x0, x1 is closer
            if (val > 0.0) {
                x1 = ans;
                qh = val;
            } else {
                x0 = ans;
                ql = val;
            }
            // stop if convergence
            if (br == 1)
                break;
        }

        return ans;
    }


    //! Calculate coefficients for cubic intrpolation spline
    //!
    //! @param vals length 4 vector of function values.
    //! @param grid length 4 vector of grid points.
    Eigen::VectorXd InterpolationGrid::find_coefs(
        const Eigen::VectorXd& vals, 
        const Eigen::VectorXd& grid
    )
    {
        Eigen::VectorXd a(4);

        double dt0 = grid(1) - grid(0);
        double dt1 = grid(2) - grid(1);
        double dt2 = grid(3) - grid(2);

        /* check for repeated points (important for boundaries) */
        if (dt1 < 1e-4) dt1 = 1.0;
        if (dt0 < 1e-4) dt0 = dt1;
        if (dt2 < 1e-4) dt2 = dt1;

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
        a(2) = -3*vals(1) + 3*vals(2) - 2*dx1 - dx2;
        a(3) = 2*vals(1) - 2*vals(2) + dx1 + dx2;

        return a;
    }

    //! Interpolate on 4 points
    //!
    //! @param x evaluation point.
    //! @param vals length 4 vector of function values.
    //! @param grid length 4 vector of grid points.
    double InterpolationGrid::interp_on_grid(
        const double& x, 
        const Eigen::VectorXd& vals, 
        const Eigen::VectorXd& grid
    )
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
    double InterpolationGrid::int_on_grid(
        const double& upr,
        const Eigen::VectorXd& vals, 
        const Eigen::VectorXd& grid
    )
    {
        int m = grid.size();
        Eigen::VectorXd tmpvals(4), tmpgrid(4), tmpa(4), a(4);
        double uprnew, newint;

        double tmpint = 0.0;

        if (upr > grid(0)) {
            // go up the grid and integrate
            for (int k = 0; k < m-1; ++k) {
                // stop loop if fully integrated
                if (upr < grid(k)) break;

                // select length 4 subvectors and calculate spline coefficients
                tmpvals(0) = vals(std::max(k-1, 0));
                tmpvals(1) = vals(k);
                tmpvals(2) = vals(k+1);
                tmpvals(3) = vals(std::min(k+2, m-1));

                tmpgrid(0) = grid(std::max(k-1, 0));
                tmpgrid(1) = grid(k);
                tmpgrid(2) = grid(k+1);
                tmpgrid(3) = grid(std::min(k+2, m-1));

                tmpa = find_coefs(tmpvals, tmpgrid);

                // don't integrate over full cell if upr is in interior
                uprnew = (upr - grid(k)) / (grid(k+1) - grid(k));
                newint = cubic_integral(0.0, fmin(1.0, uprnew), tmpa);
                tmpint += newint * (grid(k+1) - grid(k));
            }
        }

        return tmpint;
    }

    //! Inverse of the integral of a spline interpolant
    //!
    //! @param qq argument of the inverse integral (the 'quantile').
    //! @param vals vector of values to be interpolated and integrated.
    //! @param grid vector of grid points on which vals has been computed.
    //!
    //! @return Integral of interpolation spline defined by (vals, grid).
    double InterpolationGrid::inv_int_on_grid(
        const double& qq, 
        const Eigen::VectorXd& vals,
        const Eigen::VectorXd& grid
    )
    {
        int m = grid.size();
        Eigen::VectorXd tmpvals(4), tmpgrid(4), tmpa(4), a(4);
        double uprnew, newint, out, qtest;
        double tmpint = 0.0;
        int tmpk = 0;
        double q = qq;

        q *= int_on_grid(1.0, vals, grid);

        double dx = (vals(1) - vals(0)) / (grid(1) - grid(0));
        dx *= (grid(1) - grid(0));
        tmpint += vals(0)*grid(0) + (dx/2.0 * pow(grid(0), 2.0))*grid(0);
        qtest = tmpint;
        uprnew = (q - qtest)/(grid(2) - grid(1));

        // go up the grid and integrate as long as target value is above integral value
        if (q > qtest) {
            for (int k = 1; k < m-2; ++k) {
                // select length 4 subvectors and calculate spline coefficients
                tmpvals(0) = vals(k-1);
                tmpvals(1) = vals(k);
                tmpvals(2) = vals(k+1);
                tmpvals(3) = vals(k+2);

                tmpgrid(0) = grid(k-1);
                tmpgrid(1) = grid(k);
                tmpgrid(2) = grid(k+1);
                tmpgrid(3) = grid(k+2);

                tmpa = find_coefs(tmpvals, tmpgrid);
                newint = cubic_integral(0.0, 1.0, tmpa);
                tmpint += newint * (grid(k+1) - grid(k));
                tmpk = k;
                if (tmpint > q) break;
                uprnew = (q - tmpint)/(grid(tmpk+1) - grid(tmpk));
            }
        } else {
            return 2.0/dx*sqrt(q) * grid(0);
        }

        // solve in cell
        double lastgrid;
        if (tmpint > q) {
            lastgrid = grid(tmpk);
            out = lastgrid + inv_cubic_integral(uprnew, tmpa)*(grid(tmpk+1) - grid(tmpk));
        } else {
            double dx = (vals(m-1) - vals(m-2)) / (grid(m-1) - grid(m-2));
            dx *= (grid(m-1) - grid(m-2));
            out = grid(tmpk+1) + 2.0/dx*sqrt(uprnew)*(1.0 - grid(m-1));
        }
        return out;
    }
}
