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
#include <boost/python.hpp>
//#include <boost/python/numpy.hpp>

#include <bicop.hpp>

#include "git_revision.hpp"

namespace 
{
    struct bicop_wrap : Bicop, boost::python::wrapper<Bicop>
    {
        void fit(const Eigen::MatrixXd &data, std::string method) 
        { this->get_override("fit")(data, method); }

        double calculate_npars() 
        { return this->get_override("calculate_npars")(); }

        double parameters_to_tau(const Eigen::VectorXd& parameters) 
        { return this->get_override("parametrs_to_tau")(parameters); }

        Eigen::VectorXd pdf_default(const Eigen::MatrixXd& u)    
        { return this->get_override("pdf_default")(u); }

        Eigen::VectorXd hfunc1_default(const Eigen::MatrixXd& u) 
        { return this->get_override("hfunc1_default")(u); }

        Eigen::VectorXd hfunc2_default(const Eigen::MatrixXd& u) 
        { return this->get_override("hfunc2_default")(u); }

        Eigen::VectorXd hinv1_default(const Eigen::MatrixXd& u)  
        { return this->get_override("hinv1_default")(u); }

        Eigen::VectorXd hinv2_default(const Eigen::MatrixXd& u)  
        { return this->get_override("hinv2_default")(u); }
    };

//    boost::python::numpy::ndarray hello(const boost::python::numpy::ndarray &arr)
//    {
//        return arr;
//    }
}

BOOST_PYTHON_MODULE(pyvinecopulib)
{
//    Py_Initialize();
//    boost::python::numpy::initialize();

    export_git_revision();

    boost::python::class_<bicop_wrap, boost::noncopyable>("bicop", boost::python::no_init)
        .add_property("rotation", &Bicop::get_rotation, &Bicop::set_rotation)
        .add_property("parameters", &Bicop::get_parameters, &Bicop::set_parameters)
        .add_property("family", &Bicop::get_family)
    ;

    // Numpy hellow world
//    boost::python::def("hello", hello);
}

