find_package(PythonInterp REQUIRED)
add_test(NAME python_test_git_revision COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_git_revision.py)
if (WIN32)
  # https://cmake.org/pipermail/cmake/2010-December/041176.html
  # https://cmake.org/pipermail/cmake-developers/2013-June/019217.html
  find_package(Boost 1.63 COMPONENTS python REQUIRED) #TODO: reuse?
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
  get_filename_component(TMP_PYTHON_LIB_PATH ${PYTHON_LIBRARIES} DIRECTORY)

  set(PATH_STRING "$ENV{PATH};${Boost_LIBRARY_DIRS};${CMAKE_BINARY_DIR}/src/${CMAKE_BUILD_TYPE};${NLOPT_ROOT_DIR};${TMP_PYTHON_LIB_PATH}")
  string(REPLACE "\\;" ";" PATH_STRING "${PATH_STRING}")
  string(REPLACE ";" "\\;" PATH_STRING "${PATH_STRING}")

  set_property(TEST python_test_git_revision PROPERTY ENVIRONMENT "PYTHONVERBOSE=1")
  set_property(TEST python_test_git_revision APPEND PROPERTY ENVIRONMENT "PATH=${PATH_STRING}")
  set_property(TEST python_test_git_revision APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python/${CMAKE_BUILD_TYPE}")
else()
  set_tests_properties(python_test_git_revision PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python${TMP_SFX}")
endif()
