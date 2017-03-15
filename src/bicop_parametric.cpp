// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop_parametric.hpp"
#include "tools_optimization.hpp"
#include "tools_stl.hpp"
#include "tools_stats.hpp"

namespace vinecopulib
{
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
    void ParBicop::fit(const Eigen::MatrixXd &data, std::string method)
    {
        if (family_ != 0)
        {
            using namespace tools_optimization;

            std::vector<std::string> methods = {"itau", "mle"};
            if (!tools_stl::is_member(method, methods))
            {
                throw std::runtime_error("Method not implemented.");
            }

            int npars = (int) calculate_npars();
            if (method == "itau")
            {
                npars = npars - 1;
                if ((npars > 0) & (family_ != 2))
                {
                    throw std::runtime_error("itau method is not available for this family.");
                }
            }

            Eigen::MatrixXd temp_data = data;
            double tau = pairwise_ktau(temp_data);
            Eigen::VectorXd newpar = get_start_parameters(tau);

            std::string association_direction = get_association_direction();
            if (((tau < 0) & (association_direction == "positive")) | ((tau > 0) & (association_direction == "negative")))
            {
                throw std::runtime_error("The data and copula are not compatible.");
            }

            if (npars > 0)
            {
                // Create optimizer
                Optimizer optimizer(npars);

                // Set bounds and starting values
                Eigen::MatrixXd bounds = get_parameters_bounds();
                Eigen::VectorXd initial_parameters = newpar;
                ParBicopOptData my_data = {data, this, newpar(0), 0};
                if (method == "itau")
                {
                    bounds = get_parameters_bounds().block(1,0,npars,2);
                    initial_parameters = newpar.block(1,0,npars,1);
                    optimizer.set_objective(pmle_objective, &my_data);
                }
                else
                {
                    optimizer.set_objective(mle_objective, &my_data);
                }

                optimizer.set_bounds(bounds);
                Eigen::VectorXd optimized_parameters = optimizer.optimize(initial_parameters);

                if (method == "itau")
                {
                    newpar.block(1,0,npars,1) = optimized_parameters;
                }
                else
                {
                    newpar = optimized_parameters;
                }
            }
            // set the new parameters
            set_parameters(newpar);
        }

    }
}

/*void remove_row(Eigen::MatrixXd& matrix, unsigned int to_remove)
{
    unsigned int n = matrix.rows()-1;
    unsigned int m = matrix.cols();

    if(to_remove < numRows )
        matrix.block(to_remove,0,numRows-to_remove,numCols) = matrix.block(to_remove+1,0,numRows-to_remove,numCols);

    matrix.conservativeResize(numRows,numCols);
}*/
