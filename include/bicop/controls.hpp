// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <vector>
#include "family.hpp"

namespace vinecopulib {
    //! @brief A class for controlling fit of bivariate copula models.
    //!
    class ControlsBicop
    {
    public:
        // Constructor
        ControlsBicop(std::vector<BicopFamily> family_set = bicop_families::all,
                     std::string parametric_method = "mle",
                     std::string selection_criterion = "bic",
                     bool preselect_families = true);

        // Getters
        std::vector<BicopFamily> get_family_set();
        std::string get_parametric_method();
        std::string get_selection_criterion();
        bool get_preselect_families();

        // Setters
        void set_family_set(std::vector<BicopFamily> family_set);
        void set_parametric_method(std::string parametric_method);
        void set_selection_criterion(std::string selection_criterion);
        void set_preselect_families(bool preselect_families);

    private:
        std::vector<BicopFamily> family_set_;
        std::string parametric_method_;
        std::string selection_criterion_;
        bool preselect_families_;

        void check_parametric_method(std::string parametric_method);
        void check_selection_criterion(std::string selection_criterion);
    };
}
