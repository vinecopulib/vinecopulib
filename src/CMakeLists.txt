file(GLOB_RECURSE vinecopulib_sources *.cpp *.cc *c)

file(GLOB_RECURSE vinecopulib_headers ../include/*.h ../include/*.hpp)

include_directories(${external_includes} ../include)

if (BUILD_SHARED_LIBS)
    add_library(vinecopulib SHARED ${vinecopulib_sources})
else()
    add_library(vinecopulib STATIC ${vinecopulib_sources})
endif()

# TODO: remove header that should not be public
set_target_properties(vinecopulib PROPERTIES PUBLIC_HEADER "${vinecopulib_headers}")

target_link_libraries(vinecopulib ${external_libs})

set_property(TARGET vinecopulib PROPERTY POSITION_INDEPENDENT_CODE ON)

if (NOT WIN32)
    if (BUILD_SHARED_LIBS)
        install(TARGETS vinecopulib
                EXPORT vinecopulibTargets
                LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
                PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/vinecopulib")
    else()
        install(TARGETS vinecopulib
                EXPORT vinecopulibTargets
                ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
                PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/vinecopulib")
    endif()
endif()


