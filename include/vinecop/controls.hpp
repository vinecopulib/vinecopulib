// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <limits>

#include "bicop/controls.hpp"

namespace vinecopulib {
    //! @brief A class for controlling fit of bivariate copula models.
    //!
    class ControlsVinecop : public ControlsBicop
    {
    public:
        // Constructor
        ControlsVinecop();
        ControlsVinecop(std::vector<BicopFamily> family_set,
                        std::string parametric_method = "mle",
                        int truncation_level = std::numeric_limits<int>::max(),
                        std::string tree_criterion = "tau",
                        double threshold = 0.0,
                        std::string selection_criterion = "bic",
                        bool preselect_families = true,
                        bool show_trace = false);
        ControlsVinecop(const ControlsBicop controls,
                        int truncation_level = std::numeric_limits<int>::max(),
                        std::string tree_criterion = "tau",
                        double threshold = 0.0,
                        bool show_trace = false);

        // Getters
        int get_truncation_level();
        std::string get_tree_criterion();
        double get_threshold();
        bool get_show_trace();
        ControlsBicop get_controlsbicop();

        // Setters
        void set_truncation_level(int truncation_level);
        void set_tree_criterion(std::string tree_criterion);
        void set_threshold(double threshold);
        void set_show_trace(bool show_trace);
        void set_controlsbicop(ControlsBicop controls);

    private:
        int truncation_level_;
        std::string tree_criterion_;
        double threshold_;
        bool show_trace_;

        void check_truncation_level(int truncation_level);
        void check_tree_criterion(std::string tree_criterion);
        void check_threshold(double threshold);
    };
}
