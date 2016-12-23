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

#include "include/bicop_parametric.hpp"

// fit copula parameters
void ParBicop::fit(const MatXd& data, char method[])
{
    if (strcmp(method,"itau") == 0) {
        int n = data.rows();
        int d = 2;
        double tau = 0.0;
        MatXd newdata = data;
        ktau_matrix(newdata.data(), &d, &n, &tau);
        VecXd newpar = tau_to_par(tau);
        if (calculate_npars() == 1)
            set_parameters(newpar);
    } else if (strcmp(method, "mle") == 0) {

    } else {
        std::cout << "method not implemented! \n";
    }
}

// calculate number of parameters
int ParBicop::calculate_npars()
{
    // indepence copula has no parameters
    if (family_ == 0)
        return 0;
    // otherwise, return length of parameter vector
    return parameters_.size();
}
