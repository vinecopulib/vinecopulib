// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_integration.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>

namespace vinecopulib {

namespace tools_integration {

    double integrate_zero_to_one(std::function<double(double)> f)
    {
        boost::numeric::odeint::runge_kutta_dopri5<double> stepper;
        double lb = 1e-12;
        double ub = 1.0 - lb;
        double x = 0.0;
        auto ifunc = [f](const double /* x */, double &dxdt, const double t) {
            dxdt = f(t);
        };
        integrate_adaptive(boost::numeric::odeint::make_controlled(lb,
                                                                   lb,
                                                                   stepper),
                           ifunc, x, lb, ub, lb);
        return x;
    }

    Eigen::Matrix<double, Eigen::Dynamic, 2> legendre_evaluate(
            Eigen::VectorXd x) {

        Eigen::VectorXd vsub1 = x;
        Eigen::VectorXd vsub2 = Eigen::VectorXd::Ones(x.size());
        Eigen::VectorXd f = vsub2.cwiseQuotient(vsub1.cwiseAbs2() - vsub2);

        Eigen::MatrixXd vd(x.size(), 2);
        for (size_t i = 2; i <= (size_t) x.size(); ++i) {
            double ii = (double) i;
            vd.col(0) = ((2.0 * ii - 1.0) * x.cwiseProduct(vsub1) -
                         (ii - 1.0) * vsub2) / ii;
            vd.col(1) = ii * f.cwiseProduct(x.cwiseProduct(vd.col(0)) - vsub1);

            vsub2 = vsub1;
            vsub1 = vd.col(0);
        }

        return vd;
    }

    //! Get nodes and weights for the Gauss-Legendre quadrature rule
    //!
    //! @param N number of quadrature nodes and weights
    //!
    //! @return An \f$ n \times 2 \f$ matrix of with the quadrature nodes
    //! in the first column and the weights in the second
    Eigen::Matrix<double, Eigen::Dynamic, 2> legendre_rule(size_t N) {

        Eigen::VectorXd x(N);
        double NN = M_PI / ((double) N + 0.5);
        for (size_t i = 0; i < N; ++i) {
            x(i) = std::cos(NN * ((double) i + 0.75));
        }

        // value, derivative and direction
        Eigen::MatrixXd val_der = legendre_evaluate(x);
        Eigen::VectorXd dir = Eigen::VectorXd::Ones(N);
        do {
            dir = val_der.col(0).cwiseQuotient(val_der.col(1));
            x = x - dir;
            val_der = legendre_evaluate(x);
        } while (dir.cwiseAbs().maxCoeff() > 2e-16);

        Eigen::MatrixXd nodes_weights(N, 2);
        Eigen::VectorXd w = Eigen::VectorXd::Ones(N);
        nodes_weights.col(0) = (x+w)/2;
        w = (w - x.cwiseAbs2()).cwiseProduct(val_der.col(1).cwiseAbs2());
        nodes_weights.col(1) = w.cwiseInverse();

        return nodes_weights;
    }
}

}
