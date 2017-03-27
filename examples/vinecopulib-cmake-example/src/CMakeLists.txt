# Include header files
include_directories(${external_includes})

# Add vinecopulib_main executable
add_executable(vinecopulib_main vinecopulib_main.cpp)

# Link to vinecopulib
target_link_libraries(vinecopulib_main ${external_libs})