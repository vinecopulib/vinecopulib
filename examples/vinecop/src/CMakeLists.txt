# Include header files
include_directories(${external_includes})

# Add main executable
add_executable(main main.cpp)

# Link to vinecopulib and garch
target_link_libraries(main ${external_libs})