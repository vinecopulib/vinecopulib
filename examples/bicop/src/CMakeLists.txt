# Include header files
include_directories(${external_includes})

# Add vinecopulib_main executable
add_executable(main main.cpp)

# Link to vinecopulib
target_link_libraries(main ${external_libs})