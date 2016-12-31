file(GLOB_RECURSE vinecopulib_sources *.cpp *.cc *c)
file(GLOB_RECURSE vinecopulib_headers *.h *.hpp)

include_directories(${external_includes} ../include)

add_library(vinecopulib STATIC
        ${vinecopulib_sources}
        ${vinecopulib_headers})

target_link_libraries(vinecopulib ${external_libs})
