#include <boost/python.hpp>
#include <Eigen/Dense>

// https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/

namespace
{
    // TODO: change into NumPy array
    struct VecXd_to_list
    {
        static PyObject* convert(Eigen::VectorXd const &vec)
        {
            boost::python::list list;
            for (auto i=0; i < vec.size(); ++i)
                list.append(vec[i]);
            return boost::python::incref(list.ptr());
        }
    };

    struct list_to_VecXd
    {
        static void* convertible(PyObject* obj_ptr)
        {
            if (!PyList_Check(obj_ptr)) return 0;
            return obj_ptr;
        }

        static void construct(
            PyObject* obj_ptr,
            boost::python::converter::rvalue_from_python_stage1_data* data
        ) {
            void* storage = (
                (boost::python::converter::rvalue_from_python_storage<Eigen::VectorXd>*)
                data
            )->storage.bytes;

            new (storage) Eigen::VectorXd(PyList_Size(obj_ptr));

            for (Py_ssize_t i = 0; i < PyList_Size(obj_ptr); ++i)
                (*(Eigen::VectorXd*)storage)[i] = boost::python::extract<double>(PyList_GetItem(obj_ptr, i));

            data->convertible = storage;
        }
    };
}

void register_converters()
{
    boost::python::to_python_converter<Eigen::VectorXd, VecXd_to_list>();

    boost::python::converter::registry::push_back(
        &list_to_VecXd::convertible,
        &list_to_VecXd::construct,
        boost::python::type_id<Eigen::VectorXd>()
    );
}
