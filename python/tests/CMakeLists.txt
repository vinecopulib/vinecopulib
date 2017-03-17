find_package(PythonInterp REQUIRED)
foreach (test_name test_git_revision;test_bicop_class;test_vinecop_class)
  add_test(NAME ${test_name} WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.py)

  if (WIN32)
    find_package(Boost 1.63 COMPONENTS python REQUIRED) #TODO: reuse?
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
    get_filename_component(TMP_PYTHON_LIB_PATH ${PYTHON_LIBRARIES} DIRECTORY)

    # https://cmake.org/pipermail/cmake/2010-December/041176.html
    set(PATH_STRING "$ENV{PATH};${Boost_LIBRARY_DIRS};${CMAKE_BINARY_DIR}/Debug;${NLOPT_ROOT_DIR};${TMP_PYTHON_LIB_PATH}") # TODO hardcoded
    string(REPLACE "\\;" ";" PATH_STRING "${PATH_STRING}")
    string(REPLACE ";" "\\;" PATH_STRING "${PATH_STRING}")

    # https://cmake.org/pipermail/cmake-developers/2013-June/019217.html
    set_property(TEST ${test_name} PROPERTY ENVIRONMENT "PATH=${PATH_STRING}")
    set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python/src//Debug") #TODO: Debug hardcoded!
  else()
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}")
  endif()
endforeach()
