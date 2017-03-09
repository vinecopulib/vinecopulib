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

#ifndef VINECOPULIB_OPTIMIZATION_TOOLS_HPP
#define VINECOPULIB_OPTIMIZATION_TOOLS_HPP

#include <nlopt.hpp>
#include "bicop_parametric.hpp"

namespace optimization_tools {
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

    //! A class for the controls to nlopt
    //!
    class NLoptControls
    {
    public:

        //! Create controls using the default contructor
        //!
        //! @return The default NLopt controls.
        NLoptControls();

        //! Create controls by passing the arguments
        //!
        //! @param xtol_rel Relative tolerance on the parameter value
        //! @param xtol_abs Absolute tolerance on the parameter value
        //! @param ftol_rel Relative tolerance on the function value
        //! @param ftol_abs Absolue tolerance on the function value
        //! @param maxeval Maximal number of evaluations of the objective
        //! @return Custom NLopt controls.
        NLoptControls(double xtol_rel, double xtol_abs, double ftol_rel, double ftol_abs, int maxeval);


        void set_controls(nlopt::opt* opt);

        //! Getters and setters.
        //! @{
        double get_xtol_rel();
        double get_xtol_abs();
        double get_ftol_rel();
        double get_ftol_abs();
        double get_maxeval();
        //! @}
    private:
        double xtol_rel_; //! Relative tolerance on the parameter value
        double xtol_abs_; //! Absolute tolerance on the parameter value
        double ftol_rel_; //! Relative tolerance on the function value
        double ftol_abs_; //! Absolute tolerance on the function value
        int maxeval_; //! Maximal number of evaluations of the objective

        //! Sanity checks
        //! @{
        void check_parameters(double xtol_rel, double xtol_abs, double ftol_rel, double ftol_abs, int maxeval);
        //! @}
    };

    // the objective function for maximum likelihood estimation
    double mle_objective(const std::vector<double>& x,
                         std::vector<double> &,
                         void* data);

    // the objective function for profile maximum likelihood estimation
    double pmle_objective(const std::vector<double>& x,
                          std::vector<double> &,
                          void* data);

    // optimize the likelihood or profile likelihood
    std::vector<double> optimize(std::vector<double> x, nlopt::opt opt);
}
#endif //VINECOPULIB_OPTIMIZATION_TOOLS_HPP
