#include <boost/python.hpp>

#include <bicop_class.hpp>
#include <vinecop_class.hpp>

namespace 
{
    namespace bp = boost::python;

    struct bicop_wrap
    {
        BicopPtr bicop;

        bicop_wrap(const int& family, const VecXd& parameters, const int& rotation)
        {
            bicop = Bicop::create(family, parameters, rotation);
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

    std::vector<std::vector<BicopPtr>> object_to_vector_of_vector_of_BicopPtr(const bp::list& pair_copulas)
    {
        std::vector<std::vector<BicopPtr>> ret(bp::len(pair_copulas));

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

    struct vinecop_wrap : Vinecop, bp::wrapper<Vinecop>
    {
        vinecop_wrap() : Vinecop() {}

        vinecop_wrap(int d) : Vinecop(d) {}

        vinecop_wrap(const bp::list& pair_copulas, const MatXi& matrix)
            : Vinecop(
                object_to_vector_of_vector_of_BicopPtr(pair_copulas), 
                matrix
            )
        { }
    };
}

void export_bicop_class()
{
    boost::python::class_<bicop_wrap, boost::noncopyable>("bicop", boost::python::no_init)
        .def(boost::python::init<int, int>())
        .def(boost::python::init<int, VecXd, int>())
        .add_property("rotation", &bicop_wrap::get_rotation)
        .add_property("parameters", &bicop_wrap::get_parameters)
        .add_property("family", &bicop_wrap::get_family)
    ;

    //bp::to_python_converter<BicopPtr, BicopPtr_to_bicop>();
}

void export_vinecop_class()
{
    // see http://www.boost.org/doc/libs/1_63_0/libs/python/doc/html/tutorial/tutorial/functions.html#tutorial.functions.overloading
    MatXd (Vinecop::*simulate_1arg)(int) = &Vinecop::simulate;
    MatXd (Vinecop::*simulate_2arg)(int, const MatXd&) = &Vinecop::simulate;

    boost::python::class_<vinecop_wrap, boost::noncopyable>("vinecop")
        // ctors
        .def(boost::python::init<int>())
        .def(boost::python::init<boost::python::list, MatXi>()) 
        // member functions
        .def("rotation",    &Vinecop::get_rotation)
        .add_property("rotations", &Vinecop::get_rotations)
        .def("parameters",  &Vinecop::get_parameters)
        .def("family",      &Vinecop::get_family)
        .add_property("families", &Vinecop::get_families)
        .add_property("matrix",      &Vinecop::get_matrix)
        //.def("pair_copula", &Vinecop::get_pair_copula) - TODO!
        .def("pdf",         &Vinecop::pdf)
        .def("simulate", simulate_1arg) 
        .def("simulate", simulate_2arg)
        // static methods
        .def("select", &Vinecop::select, (
            boost::python::arg("data"),
            boost::python::arg("family_set"),
            boost::python::arg("method"),
            boost::python::arg("truncation_level"),
            boost::python::arg("matrix"),
            boost::python::arg("selection_criterion"),
            boost::python::arg("preselect_families"),
            boost::python::arg("show_trace")
        ))
        //.staticmethod("structure_select")
// TODO is make_pc_store needed???
//        .def("make_pair_copula_store", &Vinecop::make_pc_store)
//        .staticmethod("make_pc_store")
    ;
}
