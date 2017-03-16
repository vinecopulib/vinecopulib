#include <boost/python.hpp>

void register_eigen_converters();
void export_git_revision();
void export_bicop_class();
void export_vinecop_class();

BOOST_PYTHON_MODULE(pyvinecopulib)
{
    register_eigen_converters();
    export_git_revision();
    export_bicop_class();
    export_vinecop_class();
}
