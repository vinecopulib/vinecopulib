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

namespace optimization_tools {

    NLoptControls::NLoptControls()
    {
        xtol_rel_ = 1e-6;
        xtol_abs_ = 1e-6;
        ftol_rel_ = 1e-6;
        ftol_abs_ = 1e-6;
        maxeval_ = 1e3;
    }

    NLoptControls::NLoptControls(double xtol_rel, double xtol_abs, double ftol_rel, double ftol_abs, int maxeval)
    {
        check_parameters(xtol_rel, xtol_abs, ftol_rel, ftol_abs, maxeval);
        xtol_rel_ = xtol_rel;
        xtol_abs_ = xtol_abs;
        ftol_rel_ = ftol_rel;
        ftol_abs_ = ftol_abs;
        maxeval_ = maxeval;
    }

    void NLoptControls::check_parameters(double xtol_rel, double xtol_abs, double ftol_rel, double ftol_abs, int maxeval)
    {
        if (xtol_rel <= 0 || xtol_rel > 1)
        {
            throw std::runtime_error("xtol_rel should be in (0,1]");
        }
        if (ftol_rel <= 0 || ftol_rel > 1)
        {
            throw std::runtime_error("ftol_rel should be in (0,1]");
        }
        if (xtol_abs <= 0)
        {
            throw std::runtime_error("xtol_abs should be larger than 0");
        }
        if (ftol_abs <= 0)
        {
            throw std::runtime_error("ftol_abs should be larger than 0");
        }
        if (maxeval <= 0)
        {
            throw std::runtime_error("maxeval should be larger than 0");
        }
    }

    double NLoptControls::get_xtol_rel() {return xtol_rel_;};
    double NLoptControls::get_xtol_abs() {return xtol_abs_;};
    double NLoptControls::get_ftol_rel() {return ftol_rel_;};
    double NLoptControls::get_ftol_abs() {return ftol_abs_;};
    double NLoptControls::get_maxeval() {return maxeval_;};

    void NLoptControls::set_controls(nlopt::opt* opt)
    {
        opt->set_xtol_rel(xtol_rel_);
        opt->set_xtol_abs(xtol_abs_);
        opt->set_ftol_rel(ftol_rel_);
        opt->set_ftol_abs(ftol_abs_);
        opt->set_maxeval(maxeval_);
    };

    // the objective function for maximum likelihood estimation
    double mle_objective(const std::vector<double>& x,
                         std::vector<double> &,
                         void* data)
    {
        ParBicopMLEData* newdata = (ParBicopMLEData*) data;
        ++newdata->mle_objective_calls;
        Eigen::Map<const Eigen::VectorXd> par(&x[0], x.size());
        newdata->bicop->set_parameters(par);
        double nll = newdata->bicop->loglik(newdata->U);
        nll *= -1;
        return nll;
    }

    // the objective function for profile maximum likelihood estimation
    double pmle_objective(const std::vector<double>& x,
                          std::vector<double> &,
                          void* data)
    {
        ParBicopPMLEData* newdata = (ParBicopPMLEData*) data;
        ++newdata->pmle_objective_calls;
        VecXd par = VecXd::Ones(x.size()+1);
        par(0) = newdata->par0;
        for (unsigned int i = 0; i < x.size(); ++i)
            par(i + 1) = x[i];
        newdata->bicop->set_parameters(par);
        double nll = newdata->bicop->loglik(newdata->U);
        nll *= -1;
        return nll;
    }

    // optimize the likelihood or profile likelihood
    std::vector<double> optimize(std::vector<double> x, nlopt::opt opt)
    {
        double nll;
        std::vector<double> x_new = x;
        // optimize function
        try
        {
            opt.optimize(x_new, nll);
        } catch (nlopt::roundoff_limited err)
        {
            throw std::string("Halted because roundoff errors limited progress! ") + err.what();
        } catch (nlopt::forced_stop err)
        {
            throw std::string("Halted because of a forced termination! ") + err.what();
        } catch (std::invalid_argument err )
        {
            throw std::string("Invalid arguments. ") + err.what();
        } catch (std::bad_alloc err)
        {
            throw std::string("Ran out of memory. ") + err.what();
        } catch (std::runtime_error err)
        {
            throw std::string("Generic failure. ") + err.what();
        }
        return x_new;
    }
}
