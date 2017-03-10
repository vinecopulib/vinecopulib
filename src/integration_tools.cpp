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

#include "integration_tools.hpp"

namespace integration_tools {

    double integrate_zero_to_one(std::function<double(double)> f)
    {
        boost::numeric::odeint::runge_kutta_dopri5<double> stepper;
        double lb = 1e-12;
        double ub = 1.0 - lb;
        double x = 0.0;
        integrate_adaptive(boost::numeric::odeint::make_controlled(lb,lb,stepper),
                           [f](const double /* x */, double &dxdt, const double t) {dxdt = f(t);},
                           x, lb, ub, lb);
        return x;
    }
}
