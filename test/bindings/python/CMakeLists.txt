find_package(PythonInterp REQUIRED)
add_test(NAME python_test_git_revision COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_git_revision.py)
if (WIN32)
  # https://cmake.org/pipermail/cmake/2010-December/041176.html
  find_package(Boost 1.63 COMPONENTS python REQUIRED) #TODO: reuse?
  SET(PATH_STRING "$ENV{PATH};${Boost_LIBRARY_DIRS}")
  STRING(REPLACE "\\;" ";" PATH_STRING "${PATH_STRING}")
  STRING(REPLACE ";" "\\;" PATH_STRING "${PATH_STRING}")
  set_property(TEST python_test_git_revision PROPERTY ENVIRONMENT "PATH=${PATH_STRING}")

  set_property(TEST python_test_git_revision APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python/${CMAKE_BUILD_TYPE}")
else()
  set_tests_properties(python_test_git_revision PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python${TMP_SFX}")
endif()
