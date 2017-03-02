#include <boost/python.hpp>
//#include <boost/python/numpy.hpp>

#include <bicop.hpp>


void register_converters();
void export_git_revision();
void export_bicop_class();

BOOST_PYTHON_MODULE(pyvinecopulib)
{
//    Py_Initialize();
//    boost::python::numpy::initialize();

    register_converters();
    export_git_revision();
    export_bicop_class();
}
