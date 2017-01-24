find_package(PythonInterp REQUIRED)
add_test(NAME python_test_git_revision COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_git_revision.py)
set_tests_properties(python_test_git_revision PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python")
