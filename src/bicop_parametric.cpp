/*
    Copyright 2016 Thibault Vatter, Thomas Nagler

    This file is part of vinecopulib.

    vinecopulib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecopulib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "optimization_tools.hpp"

// calculate number of parameters
double ParBicop::calculate_npars()
{
    // indepence copula has no parameters
    if (family_ == 0)
        return 0.0;
    // otherwise, return length of parameter vector
    return (double) parameters_.size();
}

// fit
void ParBicop::fit(const MatXd &data, std::string method)
{
    if (family_ != 0)
    {
        using namespace optimization_tools;

        std::vector<std::string> methods = {"itau", "mle"};
        if (!is_member(method, methods))
        {
            throw std::runtime_error("Method not implemented.");
        }

        int npars = (int) calculate_npars();
        if (method == "itau")
        {
            npars = npars - 1;
            if (npars > 0 & family_ != 2)
            {
                throw std::runtime_error("itau method is not available for this family.");
            }
        }

        int n = data.rows();
        int d = 2;
        double tau = 0.0;
        MatXd newdata = data;
        ktau_matrix(newdata.data(), &d, &n, &tau);
        VecXd newpar = get_start_parameters(tau);

        std::string association_direction = get_association_direction();
        if (((tau < 0) & (association_direction == "positive")) | ((tau > 0) & (association_direction == "negative")))
        {
            throw std::runtime_error("The data and copula are not compatible.");
        }

        if (npars > 0)
        {
            // Create optimizer
            nlopt::opt opt(nlopt::LN_BOBYQA, npars);
            NLoptControls nlopt_controls;
            nlopt_controls.set_controls(&opt);

            // Set bounds
            std::vector<double> lb(npars);
            std::vector<double> ub(npars);
            MatXd bounds = get_parameters_bounds();
            VecXd eps = VecXd::Constant(npars,1e-6);
            if (method == "itau")
            {
                VecXd::Map(&lb[0], npars) = bounds.block(1,0,npars,1)+eps;
                VecXd::Map(&ub[0], npars) = bounds.block(1,1,npars,1)-eps;
            } else
            {
                VecXd::Map(&lb[0], npars) = bounds.col(0)+eps;
                VecXd::Map(&ub[0], npars) = bounds.col(1)-eps;
            }
            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);

            // PMLE or MLE
            if (method == "itau")
            {
                // organize data for nlopt
                MatXd U = data;
                ParBicopPMLEData my_pmle_data = {U, this, newpar(0), 0};

                // call to the optimizer
                opt.set_min_objective(pmle_objective, &my_pmle_data);

                // starting value
                std::vector<double> x(npars);
                VecXd::Map(&x[0], npars) = newpar.block(1,0,npars,1);

                // optimize function
                x = optimize(x, opt);
                for (int i = 0; i < npars; ++i)
                    newpar(i + 1) = x[i];

            }
            else
            {
                // organize data for nlopt
                MatXd U = data;
                ParBicopMLEData my_mle_data = {U, this, 0};

                // call to the optimizer
                opt.set_min_objective(mle_objective, &my_mle_data);

                // starting value
                std::vector<double> x(npars);
                VecXd::Map(&x[0], npars) = newpar;

                // optimize function
                x = optimize(x, opt);

                // save the new parameters
                Eigen::Map<const Eigen::VectorXd> parameters(&x[0], x.size());
                newpar = parameters;
            }
        }
        // set the new parameters
        set_parameters(newpar);
    }
}
