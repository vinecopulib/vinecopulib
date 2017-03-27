// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/controls.hpp"
#include "misc/tools_stl.hpp"

//! Tools for bivariate and vine copula modeling
namespace vinecopulib
{
    //! creates the controls for fitting bivariate copula models.
    //! @param family_set the set of copula families to consider (if empty, then
    //!     all families are included).
    //! @param parametric_method the fit method for parametric families;
    //!     possible choices: `"mle"`, `"itau"`.
    //! @param selection_criterion the selection criterion (`"aic"` or `"bic"`).
    //! @param preselect_families whether to exclude families before fitting
    //!     based on symmetry properties of the data.
    ControlsBicop::ControlsBicop(std::vector<BicopFamily> family_set,
                                 std::string parametric_method,
                                 std::string selection_criterion,
                                 bool preselect_families)
    {
        check_parametric_method(parametric_method);
        check_selection_criterion(selection_criterion);
        family_set_ = family_set;
        parametric_method_ = parametric_method;
        selection_criterion_ = selection_criterion;
        preselect_families_ = preselect_families;
    }

    //! Sanity checks
    //! @{
    void ControlsBicop::check_parametric_method(std::string parametric_method)
    {
        if (!tools_stl::is_member(parametric_method, {"itau", "mle"}))
        {
            throw std::runtime_error("parametric_method should be mle or itau");
        }
    }

    void ControlsBicop::check_selection_criterion(std::string selection_criterion)
    {
        if (!tools_stl::is_member(selection_criterion, {"aic", "bic"}))
        {
            throw std::runtime_error("selection_criterion should be aic or bic");
        }
    }
    //! @}

    //! Getters and setters.
    //! @{
    std::vector<BicopFamily> ControlsBicop::get_family_set()
    {
        return family_set_;
    }

    std::string ControlsBicop::get_parametric_method()
    {
        return parametric_method_;
    }

    std::string ControlsBicop::get_selection_criterion()
    {
        return selection_criterion_;
    }

    bool ControlsBicop::get_preselect_families()
    {
        return preselect_families_;
    }

    void ControlsBicop::set_family_set(std::vector<BicopFamily> family_set)
    {
        family_set_ = family_set;
    }

    void ControlsBicop::set_parametric_method(std::string parametric_method)
    {
        check_parametric_method(parametric_method);
        parametric_method_ = parametric_method;
    }

    void ControlsBicop::set_selection_criterion(std::string selection_criterion)
    {
        check_selection_criterion(selection_criterion);
        selection_criterion_ = selection_criterion;
    }

    void ControlsBicop::set_preselect_families(bool preselect_families)
    {
        preselect_families_ = preselect_families;
    }
    //! @}
}
