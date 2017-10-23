// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <stdexcept>
#include <thread>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib
{
    //! creates default controls for fitting vine copula models.
    FitControlsVinecop::FitControlsVinecop() : FitControlsBicop()
    {
        truncation_level_ = std::numeric_limits<int>::max();
        threshold_ = 0.0;
        tree_criterion_ = "tau";
        select_truncation_level_ = false;
        select_threshold_ = false;
        show_trace_ = false;
        num_threads_ = 1;
    }

    //! creates custom controls for fitting vine copula models.
    //! @param family_set see FitControlsBicop.
    //! @param parametric_method see FitControlsBicop.
    //! @param nonparametric_method see FitControlsBicop.
    //! @param nonparametric_mult see FitControlsBicop.
    //! @param truncation_level for truncated vines.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param threshold for thresholded vines (0 = no threshold).
    //! @param selection_criterion see FitControlsBicop.
    //! @param preselect_families see FitControlsBicop.
    //! @param select_truncation_level whether the truncation shall be selected
    //!     automatically.
    //! @param select_threshold whether the threshold parameter shall be 
    //!     selected automatically.
    //! @param show_trace whether to show a trace of the building progress.
    //! @param num_threads number of concurrent threads to use while fitting
    //!     pair copulas within a tree; never uses more than the number returned
    //!     by `std::thread::hardware_concurrency()``.
    FitControlsVinecop::FitControlsVinecop(std::vector<BicopFamily> family_set,
                                           std::string parametric_method,
                                           std::string nonparametric_method,
                                           double nonparametric_mult,
                                           size_t truncation_level,
                                           std::string tree_criterion,
                                           double threshold,
                                           std::string selection_criterion,
                                           bool preselect_families,
                                           bool select_truncation_level,
                                           bool select_threshold,
                                           bool show_trace,
                                           size_t num_threads) :
            FitControlsBicop(family_set,
                             parametric_method,
                             nonparametric_method,
                             nonparametric_mult,
                             selection_criterion,
                             preselect_families)
    {
        check_truncation_level(truncation_level);
        check_threshold(threshold);
        check_tree_criterion(tree_criterion);

        truncation_level_ = truncation_level;
        tree_criterion_ = tree_criterion;
        threshold_ = threshold;
        select_truncation_level_ = select_truncation_level;
        select_threshold_ = select_threshold;
        show_trace_ = show_trace;
        num_threads_ = process_num_threads(num_threads);
    }

    //! creates custom controls for fitting vine copula models.
    //! @param truncation_level for truncated vines.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param threshold for thresholded vines (0 = no threshold).
    //! @param show_trace whether to show a trace of the building progress.
    //! @param select_truncation_level whether the truncation shall be selected
    //!     automatically.
    //! @param select_threshold whether the threshold parameter shall be 
    //!     selected automatically.
    //! @param controls see FitControlsBicop.
    //! @param num_threads number of concurrent threads to use while fitting
    //!     pair copulas within a tree; never uses more than the number returned
    //!     by `std::thread::hardware_concurrency()``.
    FitControlsVinecop::FitControlsVinecop(const FitControlsBicop controls,
                                           size_t truncation_level,
                                           std::string tree_criterion,
                                           double threshold,
                                           bool select_truncation_level,
                                           bool select_threshold,
                                           bool show_trace,
                                           size_t num_threads) :
            FitControlsBicop(controls)
    {
        check_truncation_level(truncation_level);
        check_threshold(threshold);
        check_tree_criterion(tree_criterion);

        truncation_level_ = truncation_level;
        tree_criterion_ = tree_criterion;
        threshold_ = threshold;
        select_truncation_level_ = select_truncation_level;
        select_threshold_ = select_threshold;
        show_trace_ = show_trace;
        num_threads_ = process_num_threads(num_threads);
    }

    //! Sanity checks
    //! @{
    void FitControlsVinecop::check_truncation_level(size_t truncation_level)
    {
        if (truncation_level < 1) {
            throw std::runtime_error("truncation_level should greater than 1");
        }

    }
    void FitControlsVinecop::check_tree_criterion(std::string tree_criterion)
    {
        if (!tools_stl::is_member(tree_criterion, {"tau", "rho", "hoeffd"})) {
            throw std::runtime_error("tree_criterion should be tau, rho or hoeffd");
        }
    }
    void FitControlsVinecop::check_threshold(double threshold)
    {
        if (threshold < 0 || threshold > 1) {
            throw std::runtime_error("threshold should be in [0,1]");
        }
    }
    size_t FitControlsVinecop::process_num_threads(size_t num_threads)
    {
        // use at least one thread
        num_threads = std::max(num_threads, (size_t) 1);
        // don't use more threads than supported by the system
        size_t max_threads = std::thread::hardware_concurrency();
        num_threads = std::min(num_threads, max_threads);
        
        return num_threads;
    }
    //! @}

    //! Getters and setters.
    //! @{
    size_t FitControlsVinecop::get_truncation_level()
    {
        return truncation_level_;
    }

    std::string FitControlsVinecop::get_tree_criterion()
    {
        return tree_criterion_;
    }

    double FitControlsVinecop::get_threshold()
    {
        return threshold_;
    }

    bool FitControlsVinecop::get_show_trace()
    {
        return show_trace_;
    }

    bool FitControlsVinecop::get_select_truncation_level()
    {
        return select_truncation_level_;
    }
        
    bool FitControlsVinecop::get_select_threshold()
    {
        return select_threshold_;
    }
    
    size_t FitControlsVinecop::get_num_threads()
    {
        return num_threads_;
    }
    
    bool FitControlsVinecop::needs_sparse_select()
    {
        return (select_truncation_level_ | select_threshold_);
    }
    
    FitControlsBicop FitControlsVinecop::get_fit_controls_bicop()
    {
        FitControlsBicop controls_bicop(get_family_set(), 
                                        get_parametric_method(),
                                        get_nonparametric_method(),
                                        get_nonparametric_mult(),
                                        get_selection_criterion(),
                                        get_preselect_families());
        return controls_bicop;
    }

    void FitControlsVinecop::set_truncation_level(size_t truncation_level)
    {
        check_truncation_level(truncation_level);
        truncation_level_ = truncation_level;
    }

    void FitControlsVinecop::set_tree_criterion(std::string tree_criterion)
    {
        check_tree_criterion(tree_criterion);
        tree_criterion_ = tree_criterion;
    }

    void FitControlsVinecop::set_threshold(double threshold)
    {
        check_threshold(threshold);
        threshold_ = threshold;
    }

    void FitControlsVinecop::set_show_trace(bool show_trace)
    {
        show_trace_ = show_trace;
    }
    
    void FitControlsVinecop::set_select_truncation_level(bool select_truncation_level)
    {
        select_truncation_level_ = select_truncation_level;
    }
    
    void FitControlsVinecop::set_select_threshold(bool select_threshold)
    {
        select_threshold_ = select_threshold;
    }
    
    void FitControlsVinecop::set_num_threads(size_t num_threads)
    {
        num_threads_ = process_num_threads(num_threads);
    }

    void FitControlsVinecop::set_fit_controls_bicop(FitControlsBicop controls)
    {
        set_family_set(controls.get_family_set());
        set_parametric_method(controls.get_parametric_method());
        set_selection_criterion(get_selection_criterion());
        set_preselect_families(controls.get_preselect_families());
    }
    //! @}
}
