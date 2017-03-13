/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
    return rinstance_ptr->get_U();
}

VecXd FakeParBicopTest::get_parameters(RInstance *rinstance_ptr)
{
    return rinstance_ptr->get_parameters();
}

int FakeParBicopTest::get_family(RInstance *rinstance_ptr)
{
    return rinstance_ptr->get_family();
}
double FakeParBicopTest::get_tau(RInstance *rinstance_ptr)
{
    return rinstance_ptr->get_tau();
}
int FakeParBicopTest::get_rotation(RInstance *rinstance_ptr)
{
    return rinstance_ptr->get_rotation();
}
int FakeParBicopTest::get_n(RInstance *rinstance_ptr)
{
    return rinstance_ptr->get_n();
}
