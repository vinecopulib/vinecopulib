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
    double ParBicop::calculate_npars() {
        // indepence copula has no parameters
        if (family_ == BicopFamily::Indep) {
            return 0.0;
        }
        // otherwise, return length of parameter vector
        return (double) parameters_.size();
    }

    // fit
    void ParBicop::fit(const Eigen::MatrixXd &data, std::string method)
    {
        if (family_ != BicopFamily::Indep) {
            using namespace tools_optimization;

            std::vector<std::string> methods = {"itau", "mle"};
            if (!tools_stl::is_member(method, methods)) {
                throw std::runtime_error("Method not implemented.");
            }

            int npars = (int) calculate_npars();
            if (method == "itau") {
                npars = npars - 1;
                if ((npars > 0) & (family_ != BicopFamily::Student)) {
                    throw std::runtime_error("itau method is not available for this family.");
                }
            }

            auto temp_data = data;
            double tau = tools_stats::pairwise_ktau(temp_data);
            if (!tools_stl::is_member(family_, bicop_families::rotationless)) {
                if ((tau > 0) & !tools_stl::is_member(rotation_, {0, 180})) {
                    throw std::runtime_error("Copula cannot handle tau > 0");
                }
                if ((tau < 0) & !tools_stl::is_member(rotation_, {90, 270})) {
                    throw std::runtime_error("Copula cannot handle tau < 0");
                }
            }

            auto newpar = get_start_parameters(tau);
            if (npars > 0) {
                // Create optimizer
                Optimizer optimizer(npars);

                // Set bounds and starting values
                auto bounds = get_parameters_bounds();
                auto initial_parameters = newpar;
                ParBicopOptData my_data = {data, this, newpar(0), 0};
                if (method == "itau") {
                    bounds = get_parameters_bounds().block(1, 0, npars, 2);
                    initial_parameters = newpar.block(1, 0, npars, 1);
                    optimizer.set_objective(pmle_objective, &my_data);
                } else {
                    optimizer.set_objective(mle_objective, &my_data);
                }

                optimizer.set_bounds(bounds);
                auto optimized_parameters = optimizer.optimize(initial_parameters);

                if (method == "itau") {
                    newpar.block(1, 0, npars, 1) = optimized_parameters;
                } else {
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
