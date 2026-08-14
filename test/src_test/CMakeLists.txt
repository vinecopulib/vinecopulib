include_directories(${external_includes} ${CMAKE_BINARY_DIR}/generated/test/)

# Fixtures and helpers shared by several areas.
file(GLOB support_sources CONFIGURE_DEPENDS *.cpp)
list(FILTER support_sources EXCLUDE REGEX "/test_[^/]*\\.cpp$")
file(GLOB_RECURSE src_test_headers CONFIGURE_DEPENDS include/*.h include/*.hpp)

# OBJECT, not STATIC: a test's only visible effect is the static initializer
# that registers it, so a static library's members would be dropped at link
# time for want of a referenced symbol.
add_library(src_test OBJECT ${support_sources} ${src_test_headers})
target_link_libraries(src_test PUBLIC gtest vinecopulib Threads::Threads)

# One object library per area, so that editing an area recompiles only that
# area and `make test_<area>` compiles only what that area needs.
file(GLOB area_sources CONFIGURE_DEPENDS test_*.cpp)
set(test_areas "")
foreach(area_source ${area_sources})
    get_filename_component(area ${area_source} NAME_WE)
    add_library(${area} OBJECT ${area_source})
    target_link_libraries(${area} PUBLIC gtest vinecopulib Threads::Threads)
    list(APPEND test_areas ${area})
endforeach()
set(test_areas ${test_areas} PARENT_SCOPE)
