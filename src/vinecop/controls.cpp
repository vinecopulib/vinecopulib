// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "vinecop/controls.hpp"
#include "misc/tools_stl.hpp"

//! Tools for bivariate and vine copula modeling
namespace vinecopulib
{
    //! creates default controls for fitting vine copula models.
    ControlsVinecop::ControlsVinecop() : ControlsBicop()
    {
        truncation_level_ = std::numeric_limits<int>::max();
        threshold_ = 0.0;
        tree_criterion_ = "tau";
        show_trace_ = false;
    }

    //! creates custom controls for fitting vine copula models.
    //! @param family_set see ControlsBicop.
    //! @param method see ControlsBicop.
    //! @param truncation_level for truncated vines.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param threshold for thresholded vines (0 = no threshold).
    //! @param selection_criterion see ControlsBicop.
    //! @param preselect_families see ControlsBicop.
    //! @param show_trace whether to show a trace of the building progress.
    ControlsVinecop::ControlsVinecop(std::vector<BicopFamily> family_set,
                                     std::string method, int truncation_level,
                                     std::string tree_criterion,
                                     double threshold,
                                     std::string selection_criterion,
                                     bool preselect_families,
                                     bool show_trace) :
            ControlsBicop(family_set, method, selection_criterion,
                          preselect_families)
    {
        check_truncation_level(truncation_level);
        check_threshold(threshold);
        check_tree_criterion(tree_criterion);

        truncation_level_ = truncation_level;
        threshold_ = threshold;
        tree_criterion_ = tree_criterion;
        show_trace_ = show_trace;
    }

    //! creates custom controls for fitting vine copula models.
    //! @param truncation_level for truncated vines.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param threshold for thresholded vines (0 = no threshold).
    //! @param show_trace whether to show a trace of the building progress.
    //! @param controls see ControlsBicop.
    ControlsVinecop::ControlsVinecop(const ControlsBicop controls,
                                     int truncation_level,
                                     std::string tree_criterion,
                                     double threshold,
                                     bool show_trace) :
            ControlsBicop(controls)
    {
        check_truncation_level(truncation_level);
        check_threshold(threshold);
        check_tree_criterion(tree_criterion);

        truncation_level_ = truncation_level;
        threshold_ = threshold;
        tree_criterion_ = tree_criterion;
        show_trace_ = show_trace;
    }

    //! Sanity checks
    //! @{
    void ControlsVinecop::check_truncation_level(int truncation_level)
    {
        if (truncation_level < 1) {
            throw std::runtime_error("truncation_level should greater than 1");
        }

    }
    void ControlsVinecop::check_tree_criterion(std::string tree_criterion)
    {
        if (!tools_stl::is_member(tree_criterion, {"tau", "rho", "hoeffd"})) {
            throw std::runtime_error("tree_criterion should be tau, rho or hoeffd");
        }
    }
    void ControlsVinecop::check_threshold(double threshold)
    {
        if (threshold < 0 || threshold > 1) {
            throw std::runtime_error("threshold should be in [0,1]");
        }
    }
    //! @}

    //! Getters and setters.
    //! @{
    int ControlsVinecop::get_truncation_level()
    {
        return truncation_level_;
    }

    std::string ControlsVinecop::get_tree_criterion()
    {
        return tree_criterion_;
    }

    double ControlsVinecop::get_threshold()
    {
        return threshold_;
    }

    bool ControlsVinecop::get_show_trace()
    {
        return show_trace_;
    }

    ControlsBicop ControlsVinecop::get_controlsbicop()
    {
        ControlsBicop controls_bicop(get_family_set(), get_parametric_method(),
                                     get_selection_criterion(),
                                     get_preselect_families());
        return controls_bicop;
    }

    void ControlsVinecop::set_truncation_level(int truncation_level)
    {
        check_truncation_level(truncation_level);
        truncation_level_ = truncation_level;
    }

    void ControlsVinecop::set_tree_criterion(std::string tree_criterion)
    {
        check_tree_criterion(tree_criterion);
        tree_criterion_ = tree_criterion;
    }

    void ControlsVinecop::set_threshold(double threshold)
    {
        check_threshold(threshold);
        threshold_ = threshold;
    }

    void ControlsVinecop::set_show_trace(bool show_trace)
    {
        show_trace_ = show_trace;
    }

    void ControlsVinecop::set_controlsbicop(ControlsBicop controls)
    {
        set_family_set(controls.get_family_set());
        set_parametric_method(controls.get_parametric_method());
        set_selection_criterion(get_selection_criterion());
        set_preselect_families(controls.get_preselect_families());
    }
    //! @}
}
