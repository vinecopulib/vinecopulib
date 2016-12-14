//
// Created by Vatter Thibault on 13/12/16.
//

#include "include/par_bicop.h"

// getters --------------------------------

int ParBicop::get_family()
{
    return this->family_;
}

VecXd ParBicop::get_parameters()
{
    return this->parameters_;
}

// setters --------------------------------

void ParBicop::set_family(int family)
{
    family_ = family;
}

void ParBicop::set_parameters(VecXd parameters)
{
    parameters_ = parameters;
}

// Family-specific functions --------------------------------

// calculate number of parameters
int ParBicop::calculate_npars()
{
    // indepence copula has no parameters
    if (family_ == 0)
        return(0);

    // otherwise, return length of parameter vector
    return(int(parameters_.size()));
}
