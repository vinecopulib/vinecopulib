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
# NumPy support was introduced in 1.63
if (APPLE OR WIN32) # Windows & Homebrew convention
  set(TMP_BOOST_PTHON_COMPONENT "python")
else() # Debian convention
  set(TMP_BOOST_PTHON_COMPONENT "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
endif()
find_package(Boost 1.63 COMPONENTS ${TMP_BOOST_PYTHON_COMPONENT} REQUIRED)
unset(TMP_BOOST_PYTHON_COMPONENT)
#<TEMP!>
message(${Boost_LIBRARIES})
#</TEMP!>
target_link_libraries(pyvinecopulib ${Boost_LIBRARIES})
target_include_directories(pyvinecopulib PUBLIC ${Boost_INCLUDE_DIRS})


# Eigen
#TODO: find_package() ...
target_include_directories(pyvinecopulib PUBLIC ${EIGEN3_INCLUDE_DIR})


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
