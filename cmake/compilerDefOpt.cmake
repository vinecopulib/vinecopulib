# Append, never assign: assigning shadows the cache entry CMake initializes from
# $CXXFLAGS, discarding distribution hardening flags and user-supplied options.

if(NOT WIN32)

    if(VINECOPULIB_STRICT_COMPILER AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        string(APPEND CMAKE_CXX_FLAGS " -Werror -Wno-delete-non-virtual-dtor -Wall -Wextra -Wpedantic -Wconversion")
        string(APPEND CMAKE_CXX_FLAGS " -Wstrict-aliasing -pedantic -fmax-errors=5 -Werror=return-type")
        string(APPEND CMAKE_CXX_FLAGS " -Wunreachable-code -Wcast-align -Wcast-qual")
        string(APPEND CMAKE_CXX_FLAGS " -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op")
        string(APPEND CMAKE_CXX_FLAGS " -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual")
        string(APPEND CMAKE_CXX_FLAGS " -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel")
        string(APPEND CMAKE_CXX_FLAGS " -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option")
        string(APPEND CMAKE_CXX_FLAGS " -Wbool-operation")
        string(APPEND CMAKE_CXX_FLAGS " -Wbuiltin-macro-redefined -Wundef")
    else()
        string(APPEND CMAKE_CXX_FLAGS " -Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type")
    endif()

    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -g -O0 -DDEBUG")

    # Opt-in: the result only runs on CPUs at least as new as the builder's.
    if(VINECOPULIB_NATIVE_ARCH)
        if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64" AND CMAKE_SYSTEM_NAME STREQUAL "Darwin")
            string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mcpu=apple-m1")
        else()
            string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native")
        endif()
    endif()

    if(VINECOPULIB_SANITIZERS)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    endif()

    if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING AND NOT WIN32 AND VINECOPULIB_CODE_COVERAGE)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fprofile-arcs -ftest-coverage")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

        if(NOT EXISTS ${CMAKE_CXX_COMPILER})
            message( FATAL_ERROR "Clang++ not found. " )
        endif()

        string(APPEND CMAKE_CXX_FLAGS " -Wno-deprecated-register")
        string(APPEND CMAKE_CXX_FLAGS " -Wno-unused-const-variable")
        string(APPEND CMAKE_CXX_FLAGS " -fcolor-diagnostics")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 14.0)
            string(APPEND CMAKE_CXX_FLAGS " -Wbitwise-instead-of-logical")
        endif ()
    endif()
endif()

if (MSVC)
    string(APPEND CMAKE_CXX_FLAGS " /bigobj")
endif()

set(VINECOPULIB_DEFINITIONS
  BOOST_NO_AUTO_PTR
  BOOST_ALLOW_DEPRECATED_HEADERS
  BOOST_MATH_PROMOTE_DOUBLE_POLICY=false
  BOOST_ALL_NO_LIB
  # Selects Boost's mt19937 over std::mt19937 in wdm's RandomGenerator
  # (wdm/random.hpp). Changing it changes the "random" ties method in
  # to_pseudo_obs, and so every tll fit.
  USE_BOOST
)

add_compile_definitions(${VINECOPULIB_DEFINITIONS})
