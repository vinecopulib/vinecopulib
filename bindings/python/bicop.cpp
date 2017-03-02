#include <boost/python.hpp>
//#include <boost/python/numpy.hpp>

#include <bicop_class.hpp>

namespace 
{
    struct bicop_wrap
    {
        BicopPtr bicop;

        bicop_wrap(const int& family, const boost::python::object& parameters, const int& rotation)
        {
            bicop = Bicop::create(family, boost::python::extract<VecXd>(parameters)(), rotation);
        }

        bicop_wrap(const int& family, const int& rotation)
        {
            bicop = Bicop::create(family, rotation);
        }

        int get_family() 
        { 
            return bicop->get_family(); 
        }

        int get_rotation() 
        { 
            return bicop->get_rotation(); 
        }

        VecXd get_parameters() 
        {
            return bicop->get_parameters();
        }
    };
}

void export_bicop_class()
{
    boost::python::class_<bicop_wrap, boost::noncopyable>("bicop", boost::python::no_init)
        .def(boost::python::init<int, int>())
        .def(boost::python::init<int, boost::python::object, int>())
        .add_property("rotation", &bicop_wrap::get_rotation)
        .add_property("parameters", &bicop_wrap::get_parameters)
        .add_property("family", &bicop_wrap::get_family)
    ;
}
