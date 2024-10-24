# Add main executable
add_executable(main main.cpp)

# Link to vinecopulib
target_link_libraries(main ${VINECOPULIB_LIBRARIES})
