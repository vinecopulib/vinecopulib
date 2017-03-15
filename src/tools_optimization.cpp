// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "tools_optimization.hpp"

namespace tools_optimization {

    Optimizer::Optimizer(unsigned int n_parameters)
    {
        if (n_parameters < 1)
        {
            throw std::runtime_error("n_parameters should be larger than 0.");
        }
        n_parameters_ = n_parameters;
        opt_ = nlopt::opt(nlopt::LN_BOBYQA, n_parameters);
        controls_ = NLoptControls();
        controls_.set_controls(&opt_);
    }

    Optimizer::Optimizer(unsigned int n_parameters, double xtol_rel, double xtol_abs,
                         double ftol_rel, double ftol_abs, int maxeval)
    {
        if (n_parameters < 1)
        {
            throw std::runtime_error("n_parameters should be larger than 0.");
        }
        n_parameters_ = n_parameters;
        opt_ = nlopt::opt(nlopt::LN_BOBYQA, n_parameters);
        controls_ = NLoptControls(xtol_rel, xtol_abs, ftol_rel, ftol_abs, maxeval);
        controls_.set_controls(&opt_);
    }

    void Optimizer::set_bounds(vinecopulib::MatXd bounds)
    {
        if (bounds.rows() != n_parameters_ || bounds.cols() != 2)
        {
            throw std::runtime_error("Bounds should be a two column matrix with n_parameters_ rows.");
        }

        std::vector<double> lb(n_parameters_);
        std::vector<double> ub(n_parameters_);
        vinecopulib::VecXd eps = vinecopulib::VecXd::Constant(n_parameters_,1e-6);
        vinecopulib::VecXd::Map(&lb[0], n_parameters_) = bounds.col(0)+eps;
        vinecopulib::VecXd::Map(&ub[0], n_parameters_) = bounds.col(1)-eps;
        opt_.set_lower_bounds(lb);
        opt_.set_upper_bounds(ub);
    }

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
                         std::vector<double>&,
                         void* f_data)
    {
        ParBicopOptData* newdata = (ParBicopOptData*) f_data;
        ++newdata->objective_calls;
        Eigen::Map<const Eigen::VectorXd> par(&x[0], x.size());
        newdata->bicop->set_parameters(par);
        double nll = newdata->bicop->loglik(newdata->U);
        nll *= -1;
        return nll;
    }

    // the objective function for profile maximum likelihood estimation
    double pmle_objective(const std::vector<double>& x,
                          std::vector<double>&,
                          void* f_data)
    {
        ParBicopOptData* newdata = (ParBicopOptData*) f_data;
        ++newdata->objective_calls;
        vinecopulib::VecXd par = vinecopulib::VecXd::Ones(x.size()+1);
        par(0) = newdata->par0;
        for (unsigned int i = 0; i < x.size(); ++i)
            par(i + 1) = x[i];
        newdata->bicop->set_parameters(par);
        double nll = newdata->bicop->loglik(newdata->U);
        nll *= -1;
        return nll;
    }

    void Optimizer::set_objective(nlopt::vfunc f, void* f_data)
    {
        opt_.set_min_objective(f, f_data);
    }

    // optimize the likelihood or profile likelihood
    vinecopulib::VecXd Optimizer::optimize(vinecopulib::VecXd initial_parameters)
    {
        if (initial_parameters.size() != n_parameters_)
        {
            throw std::string("The size of x should be n_parameters_.");
        }

        double nll;
        std::vector<double> x(n_parameters_);
        vinecopulib::VecXd::Map(&x[0], n_parameters_) = initial_parameters;
        try
        {
            opt_.optimize(x, nll);
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
        } catch (...) {}

        Eigen::Map<const Eigen::VectorXd> optimized_parameters(&x[0], x.size());
        return optimized_parameters;
    }
}
