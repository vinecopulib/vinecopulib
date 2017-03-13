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

#include "include/parbicop_test.hpp"

void FakeParBicopTest::set_family(int family, int rotation)
{
    std::vector<int> rotation_less_fams = {0, 1, 2, 5};
    if (is_member(family, rotation_less_fams) || rotation == 0)
    {
        family_ = family;
    }
    else
    {
        if (rotation == 90)
        {
            family_ = family + 20;
        }
        else if (rotation == 180)
        {
            family_ = family + 10;
        }
        else
        {
            family_ = family + 30;
        }
    }
}
void FakeParBicopTest::set_parameters(VecXd parameters)
{
    if (parameters.size() > 0)
    {
        par_ = parameters(0);
    }
    if (parameters.size() > 1)
    {
        par2_ = parameters(1);
    }
}
void FakeParBicopTest::set_n(int n)
{
    n_ = n;
}
int FakeParBicopTest::get_family()
{
    return family_;
}
int FakeParBicopTest::get_n()
{
    return n_;
}
double FakeParBicopTest::get_par()
{
    return par_;
}
double FakeParBicopTest::get_par2()
{
    return par2_;
}