add_library(pyvinecopulib SHARED git_revision.cpp bicop.cpp)
add_dependencies(pyvinecopulib git_revision.h)
set_target_properties(pyvinecopulib PROPERTIES SUFFIX ".so") # e.g. Mac defaults to .dylib which is not looked for by Python
target_link_libraries(pyvinecopulib vinecopulib)
target_include_directories(pyvinecopulib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# informing the Python bindings where to find Python
find_package(PythonInterp REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
target_include_directories(pyvinecopulib PUBLIC ${PYTHON_INCLUDE_DIRS})
target_link_libraries(pyvinecopulib ${PYTHON_LIBRARIES})


# informing the Python bindings where to find Boost.Python 
# TODO! NumPy support was introduced in 1.63
find_package(Boost 1.54 COMPONENTS python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} QUIET REQUIRED)
message(WARNING "FIXME! You may safely ignore warning the 'No header definder for python-py...' if present! (see https://gitlab.kitware.com/cmake/cmake/issues/16335)")
target_link_libraries(pyvinecopulib ${Boost_LIBRARIES})


# installation (see http://www.cmake.org/pipermail/cmake/2013-March/054138.html)
if (NOT PYTHONLIBS_VERSION_STRING EQUAL PYTHON_VERSION_STRING)
    message(FATAL_ERROR "Python libs version (${PYTHONLIBS_VERSION_STRING}) does not match version of the Python interpreter (${PYTHON_VERSION_STRING})")
endif()
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
   	"import site, sys; sys.stdout.write(site.PREFIXES[-1])"
   OUTPUT_VARIABLE PYTHON_PREFIX
)
file ( TO_CMAKE_PATH "${PYTHON_PREFIX}" PYTHON_PREFIX )
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
   	"import site, sys; sys.stdout.write(site.getsitepackages()[-1])"
   OUTPUT_VARIABLE PYTHON_SITE_DIR
)
file ( TO_CMAKE_PATH "${PYTHON_SITE_DIR}" PYTHON_SITE_DIR )
string ( REGEX REPLACE "^${PYTHON_PREFIX}/" ""
   PYTHON_SITE_DIR "${PYTHON_SITE_DIR}"
)

install ( TARGETS pyvinecopulib
   LIBRARY
     DESTINATION ${PYTHON_SITE_DIR}
     COMPONENT library
)
