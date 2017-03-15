include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(src)
if(BUILD_TESTING)
    add_subdirectory(test)
endif(BUILD_TESTING)
#add_subdirectory(bindings)
