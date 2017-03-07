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

#ifndef VINECOPULIB_BICOP_PARAMETRIC_HPP
#define VINECOPULIB_BICOP_PARAMETRIC_HPP

#include "bicop_class.hpp"

class ParBicop : public Bicop {

public:
    // fit copula parameters
    void fit(const MatXd& data, std::string method);

    // link between Kendall's tau and the par_bicop parameter
    virtual VecXd tau_to_parameters(const double& tau) = 0;
    virtual double parameters_to_tau(const VecXd& parameters) = 0;
    double calculate_tau() {return this->parameters_to_tau(parameters_);}

    // number of parameters
    double calculate_npars();

private:
    virtual VecXd get_start_parameters(const double tau) = 0;
};


//! A helper struct for nlopt maximum likelihood estimation
//!
typedef struct
{
    MatXd& U; //! The data
    ParBicop* bicop; //! A pointer to the bivariate copula to optimize
    unsigned int mle_objective_calls; //! The number of evaluations of the objective
} ParBicopMLEData;

//! A helper struct for nlopt profile maximum likelihood estimation
//!
typedef struct
{
    MatXd& U; //! The data
    ParBicop* bicop; //! A pointer to the bivariate copula to optimize
    double par0;  //! The main dependence parameter
    unsigned int pmle_objective_calls; //! The number of evaluations of the objective
} ParBicopPMLEData;

double mle_objective(const std::vector<double>& x,
                    std::vector<double> &,
                    void* data);
double pmle_objective(const std::vector<double>& x,
                      std::vector<double> &,
                      void* data);

#endif
