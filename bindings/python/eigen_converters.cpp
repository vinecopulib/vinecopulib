#include <boost/python.hpp>
#include <numpy/ndarrayobject.h> 
#include <Eigen/Dense>

// https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
// http://stackoverflow.com/questions/10701514/how-to-return-numpy-array-from-boostpython

#include <iostream>

namespace
{
    namespace bp = boost::python;

    struct VecXd_to_list
    {
        static PyObject* convert(Eigen::VectorXd const &vec)
        {
            bp::list list;
            for (auto i=0; i < vec.size(); ++i)
                list.append(vec[i]);
            return bp::incref(list.ptr());
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
            bp::converter::rvalue_from_python_stage1_data* data
        ) {
            void* storage = (
                (bp::converter::rvalue_from_python_storage<Eigen::VectorXd>*) data
            )->storage.bytes;

            new (storage) Eigen::VectorXd(PyList_Size(obj_ptr));

            for (Py_ssize_t i = 0; i < PyList_Size(obj_ptr); ++i)
                (*(Eigen::VectorXd*)storage)[i] = bp::extract<double>(PyList_GetItem(obj_ptr, i));

            data->convertible = storage;
        }
    };

    struct MatXi_to_ndarray
    {
        static PyObject* convert(Eigen::MatrixXi const &mat)
        {
            npy_intp size[2];
            size[1] = mat.rows();
            size[0] = mat.cols();

            Eigen::MatrixXi::Scalar * data = const_cast<Eigen::MatrixXi::Scalar *>(&mat(0,0)); 

// TODO: assert sizeof(Eigen::MatrixXi::Scalar) == sizof(NPY_INT)
// TODO: assert strides OK
// TODO: handle zero size

            PyObject * pyObj = PyArray_SimpleNewFromData(2, size, NPY_INT32, data);
            bp::handle<> handle( pyObj );
            bp::numeric::array arr( handle );

            return bp::incref(arr.copy().ptr());
        }
    };

    template<int NumpyType, class EigenType>
    struct ndarray_to_MatX
    {
        static void* convertible(PyObject* obj_ptr)
        {
// TODO: check if contiguous
// TODO: should we throw here?
// TODO: assert sizeof
            if (PyArray_NDIM(obj_ptr) != 2) return 0;
            if (PyArray_TYPE(obj_ptr) != NumpyType) return 0;
            return obj_ptr;
        }

        static void construct(
            PyObject* obj_ptr,
            bp::converter::rvalue_from_python_stage1_data* data
        ) {
            void* storage = (
                (bp::converter::rvalue_from_python_storage<EigenType>*) data
            )->storage.bytes;

            new (storage) EigenType(
                PyArray_DIM(obj_ptr, 1), 
                PyArray_DIM(obj_ptr, 0)
            );

            for (Py_ssize_t i = 0; i < PyArray_DIM(obj_ptr, 0); ++i)
                for (Py_ssize_t j = 0; j < PyArray_DIM(obj_ptr, 1); ++j)
                    (*(EigenType*)storage)(j, i) = *(typename EigenType::Scalar*)(PyArray_GETPTR2(obj_ptr, i, j));

            data->convertible = storage;
        }
    };
}

void register_eigen_converters()
{
    bp::numeric::array::set_module_and_type("numpy", "ndarray");
    _import_array();

    bp::to_python_converter<Eigen::VectorXd, VecXd_to_list>();

    bp::to_python_converter<Eigen::MatrixXi, MatXi_to_ndarray>();   

    bp::converter::registry::push_back(
        &list_to_VecXd::convertible,
        &list_to_VecXd::construct,
        bp::type_id<Eigen::VectorXd>()
    );

    bp::converter::registry::push_back(
        &ndarray_to_MatX<NPY_INT32, Eigen::MatrixXi>::convertible,
        &ndarray_to_MatX<NPY_INT32, Eigen::MatrixXi>::construct,
        bp::type_id<Eigen::MatrixXi>()
    );

    bp::converter::registry::push_back(
        &ndarray_to_MatX<NPY_DOUBLE, Eigen::MatrixXd>::convertible,
        &ndarray_to_MatX<NPY_DOUBLE, Eigen::MatrixXd>::construct,
        bp::type_id<Eigen::MatrixXd>()
    );
}
