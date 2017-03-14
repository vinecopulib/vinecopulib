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

#include "bicop_parametric.hpp"
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

/*void remove_row(MatXd& matrix, unsigned int to_remove)
{
    unsigned int n = matrix.rows()-1;
    unsigned int m = matrix.cols();

    if(to_remove < numRows )
        matrix.block(to_remove,0,numRows-to_remove,numCols) = matrix.block(to_remove+1,0,numRows-to_remove,numCols);

    matrix.conservativeResize(numRows,numCols);
}*/

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
            if ((npars > 0) & (family_ != 2))
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
            Optimizer optimizer(npars);

            // Set bounds and starting values
            MatXd bounds = get_parameters_bounds();
            VecXd initial_parameters = newpar;
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
            VecXd optimized_parameters = optimizer.optimize(initial_parameters);

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
