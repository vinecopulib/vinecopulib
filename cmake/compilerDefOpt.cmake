if(NOT WIN32)

    if(STRICT_COMPILER AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
      set(CMAKE_CXX_FLAGS                "-Werror -Wno-delete-non-virtual-dtor -Wall -Wextra -Wconversion")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wstrict-aliasing -pedantic -fmax-errors=5 -Werror=return-type")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wunreachable-code -Wcast-align -Wcast-qual")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wno-variadic-macros -Wno-parentheses  -fdiagnostics-show-option")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wbool-operation")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wbuiltin-macro-redefined -Wundef")
    else()
        set(CMAKE_CXX_FLAGS                "-Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type")
    endif()

    set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0 -DDEBUG ")
    
    # Detect Apple M1 and apply specific flags
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64" AND CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        message(STATUS "Configuring for Apple M1")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -arch arm64 -mcpu=apple-m1 -DNDEBUG")
    else()
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
    endif()

    if(OPT_ASAN)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
        # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread -fno-omit-frame-pointer")
    endif()

    if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING AND NOT WIN32 AND CODE_COVERAGE)
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

        if(NOT EXISTS ${CMAKE_CXX_COMPILER})
            message( FATAL_ERROR "Clang++ not found. " )
        endif()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register") # -Qunused-arguments
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-const-variable ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 14.0)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wc++17-attribute-extensions")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wbitwise-instead-of-logical")
        endif ()
    endif()
endif()

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

add_compile_definitions(
  BOOST_NO_AUTO_PTR
  BOOST_ALLOW_DEPRECATED_HEADERS
  BOOST_MATH_PROMOTE_DOUBLE_POLICY=false
  BOOST_ALL_NO_LIB
  USE_BOOST
)