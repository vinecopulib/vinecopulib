# Include header files
include_directories(${external_includes})

# Add main executable
add_executable(main main.cpp)

# Link to vinecopulib if vinecopulib has been built as a shared lib
# and to pthreads on some linux systems (does nothing otherwise)
target_link_libraries(main ${VINECOPULIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})