#include <boost/python.hpp>

#include <bicop.hpp>

#include "git_revision.hpp"

namespace 
{
    struct bicop_wrap : Bicop, boost::python::wrapper<Bicop>
    {
        void fit(const MatXd &data, std::string method) 
        { this->get_override("fit")(data, method); }

        double calculate_npars() 
        { return this->get_override("calculate_npars")(); }

        double parameters_to_tau(const VecXd& parameters) 
        { return this->get_override("parametrs_to_tau")(parameters); }

        VecXd pdf_default(const MatXd& u)    
        { return this->get_override("pdf_default")(u); }

        VecXd hfunc1_default(const MatXd& u) 
        { return this->get_override("hfunc1_default")(u); }

        VecXd hfunc2_default(const MatXd& u) 
        { return this->get_override("hfunc2_default")(u); }

        VecXd hinv1_default(const MatXd& u)  
        { return this->get_override("hinv1_default")(u); }

        VecXd hinv2_default(const MatXd& u)  
        { return this->get_override("hinv2_default")(u); }
    };
}

BOOST_PYTHON_MODULE(libpyvinecopulib)
{
    export_git_revision();

    boost::python::class_<bicop_wrap, boost::noncopyable>("bicop", boost::python::no_init)
        .add_property("rotation", &Bicop::get_rotation, &Bicop::set_rotation)
    ;
}

