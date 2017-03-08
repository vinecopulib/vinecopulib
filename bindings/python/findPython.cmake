
find_package(PythonInterp REQUIRED)


if(APPLE)
    find_file(_find_lib_python_py findPython.py PATHS ${CMAKE_CURRENT_LIST_DIR})
    execute_process(COMMAND ${PYTHON_EXECUTABLE}  ${_find_lib_python_py} OUTPUT_VARIABLE python_config)
    string(REGEX REPLACE ".*\nshort_version:([^\n]+).*$" "\\1" PYTHON_SHORT_VERSION ${python_config})
    string(REGEX REPLACE ".*\npy_inc_dir:([^\n]+).*$" "\\1" PYTHON_INCLUDE_DIR ${python_config})
    string(REGEX REPLACE "include/python${PYTHON_SHORT_VERSION}" "Python" PYTHON_LIBRARY ${PYTHON_INCLUDE_DIR})
    if(NOT PYTHON_LIBRARY)
        message(FATAL_ERROR "Could not find Python.")
    endif()
    set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})
else(APPLE)
    find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
endif(APPLE)