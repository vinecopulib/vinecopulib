#include <boost/python.hpp>

#include <git_revision.h>

BOOST_PYTHON_MODULE(libpyvinecopulib)
{
    namespace bp = boost::python;

    // specify that this module is actually a package
    bp::object package = bp::scope();
    package.attr("__path__") = "libcloudphxx";

    // exposing git revision id
    package.attr("git_revision") = GIT_REVISION;
}
