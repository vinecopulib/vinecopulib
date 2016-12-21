/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecoplib.

    vinecoplib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecoplib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "include/bicop_parametric_test.hpp"

VecXd FakeParBicopTest::eval_in_R(RInstance *rinstance_ptr, std::string eval_fct, int start)
{
    VecXd result = rinstance_ptr->eval_in_R(eval_fct, start);
    return(result);
}

void FakeParBicopTest::set_family(RInstance *rinstance_ptr, int family)
{
    rinstance_ptr->set_family(family);
}
void FakeParBicopTest::set_rotation(RInstance *rinstance_ptr, int rotation)
{
    rinstance_ptr->set_rotation(rotation);
}
void FakeParBicopTest::set_tau(RInstance *rinstance_ptr, double tau)
{
    rinstance_ptr->set_tau(tau);
}
void FakeParBicopTest::set_parameters(RInstance *rinstance_ptr, VecXd parameters)
{
    rinstance_ptr->set_parameters(parameters);
}
void FakeParBicopTest::change_n(RInstance *rinstance_ptr, int n)
{
    rinstance_ptr->change_n(n);
}
MatXd FakeParBicopTest::get_U(RInstance *rinstance_ptr)
{
    return(rinstance_ptr->get_U());
}

VecXd FakeParBicopTest::get_parameters(RInstance *rinstance_ptr)
{
    return(rinstance_ptr->get_parameters());
}

int FakeParBicopTest::get_family(RInstance *rinstance_ptr)
{
    return(rinstance_ptr->get_family());
}
double FakeParBicopTest::get_tau(RInstance *rinstance_ptr)
{
    return(rinstance_ptr->get_tau());
}
int FakeParBicopTest::get_rotation(RInstance *rinstance_ptr)
{
    return(rinstance_ptr->get_rotation());
}