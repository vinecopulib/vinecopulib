/*
* The MIT License (MIT)
*
* Author: Sylwester Arabas, Ahmad Farhat
* Copyright © 2017 Chatham Financial
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <boost/python.hpp>

#include <bicop/class.hpp>
#include <vinecop/class.hpp>

namespace 
{
    namespace bp = boost::python;

    struct bicop_wrap
    {
        vinecopulib::BicopPtr bicop;

        bicop_wrap(const vinecopulib::BicopFamily& family, const int& rotation, const Eigen::VectorXd& parameters)
        {
            bicop = vinecopulib::Bicop::create(family, rotation, parameters);
        }

        bicop_wrap(const vinecopulib::BicopFamily& family, const int& rotation)
        {
            bicop = vinecopulib::Bicop::create(family, rotation);
        }

        bicop_wrap(const vinecopulib::BicopFamily& family)
        {
            bicop = vinecopulib::Bicop::create(family);
        }

        vinecopulib::BicopFamily get_family() 
        { 
            return bicop->get_family(); 
        }

        int get_rotation() 
        { 
            return bicop->get_rotation(); 
        }

        Eigen::VectorXd get_parameters() 
        {
            return bicop->get_parameters();
        }
    };

    std::vector<std::vector<vinecopulib::BicopPtr>> object_to_vector_of_vector_of_BicopPtr(const bp::list& pair_copulas)
    {
        std::vector<std::vector<vinecopulib::BicopPtr>> ret(bp::len(pair_copulas));

        for (auto i=0; i<ret.size(); ++i)
        {
            auto list = bp::extract<bp::list>(pair_copulas[i])();
            ret[i].resize(bp::len(list));
            for (auto j=0; j<bp::len(list); ++j)
            {
                auto elem = list[j];
                ret[i][j] = bp::extract<bicop_wrap>(elem)().bicop;
            }
        }

        return ret;
    }

    struct vinecop_wrap : vinecopulib::Vinecop, bp::wrapper<vinecopulib::Vinecop>
    {
        vinecop_wrap() : vinecopulib::Vinecop() {}

        vinecop_wrap(int d) : vinecopulib::Vinecop(d) {}

        vinecop_wrap(const bp::list& pair_copulas, const Eigen::MatrixXi& matrix)
            : vinecopulib::Vinecop(
                object_to_vector_of_vector_of_BicopPtr(pair_copulas), 
                matrix
            )
        { }
    };
}

void export_family_enums()
{
    boost::python::enum_<vinecopulib::BicopFamily>("BicopFamily")
        .value("indep",    vinecopulib::BicopFamily::indep)
        .value("gaussian", vinecopulib::BicopFamily::gaussian)
        .value("student",  vinecopulib::BicopFamily::student)
        .value("clayton",  vinecopulib::BicopFamily::clayton)
        .value("gumbel",   vinecopulib::BicopFamily::gumbel)
        .value("frank",    vinecopulib::BicopFamily::frank)
        .value("joe",      vinecopulib::BicopFamily::joe)
        .value("bb1",      vinecopulib::BicopFamily::bb1)
        .value("bb6",      vinecopulib::BicopFamily::bb6)
        .value("bb7",      vinecopulib::BicopFamily::bb7)
        .value("bb8",      vinecopulib::BicopFamily::bb8)
        .value("tll0",     vinecopulib::BicopFamily::tll0)
        .export_values()
        ;
}

void export_bicop_class()
{
    boost::python::class_<bicop_wrap, boost::noncopyable>("bicop", boost::python::no_init)
        .def(boost::python::init<vinecopulib::BicopFamily>())
        .def(boost::python::init<vinecopulib::BicopFamily, int>())
        .def(boost::python::init<vinecopulib::BicopFamily, int, Eigen::VectorXd>())
        .add_property("rotation", &bicop_wrap::get_rotation)
        .add_property("parameters", &bicop_wrap::get_parameters)
        .add_property("family", &bicop_wrap::get_family)
    ;

    //bp::to_python_converter<BicopPtr, BicopPtr_to_bicop>();
}

void export_vinecop_class()
{
    boost::python::class_<vinecop_wrap, boost::noncopyable>("vinecop")
        // ctors
        .def(boost::python::init<int>())
        .def(boost::python::init<boost::python::list, Eigen::MatrixXi>()) 
        // member functions
        .def("rotation",    &vinecopulib::Vinecop::get_rotation)
        .add_property("all_rotations", &vinecopulib::Vinecop::get_all_rotations)
        .def("parameters",  &vinecopulib::Vinecop::get_parameters)
        .def("all_parameters",  &vinecopulib::Vinecop::get_all_parameters)
        .def("family",      &vinecopulib::Vinecop::get_family)
        .add_property("all_families", &vinecopulib::Vinecop::get_all_families)
        .add_property("matrix",       &vinecopulib::Vinecop::get_matrix)
        //.def("pair_copula", &vinecopulib::Vinecop::get_pair_copula) - TODO!
        .def("pdf",         &vinecopulib::Vinecop::pdf)
        .def("simulate", &vinecopulib::Vinecop::simulate) 
        .def("inverse_rosenblatt", &vinecopulib::Vinecop::inverse_rosenblatt)
        // static methods
        .def("select", &vinecopulib::Vinecop::select, (
            boost::python::arg("data"),
            boost::python::arg("family_set"),
            boost::python::arg("method"),
            boost::python::arg("truncation_level"),
            boost::python::arg("matrix"),
            boost::python::arg("selection_criterion"),
            boost::python::arg("preselect_families"),
            boost::python::arg("show_trace")
        )) // TODO: bp::return_value_policy<bp::manage_new_object>()
        //.staticmethod("structure_select")
// TODO is make_pc_store needed???
//        .def("make_pair_copula_store", &vinecopulib::Vinecop::make_pc_store)
//        .staticmethod("make_pc_store")
    ;
}
