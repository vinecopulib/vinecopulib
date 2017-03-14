/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
        if (!stl_tools::is_member(method, methods))
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
