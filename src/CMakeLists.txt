file(GLOB_RECURSE vinecopulib_sources *.cpp *.cc *c)
file(GLOB_RECURSE vinecopulib_headers include/*.h include/*.hpp)

include_directories(${external_includes})

add_library(vinecopulib STATIC
        ${vinecopulib_sources}
        ${vinecopulib_headers})

target_link_libraries(vinecopulib ${external_libs})
