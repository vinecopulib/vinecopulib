if(NOT WIN32)

    if(STRICT_COMPILER AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS                "-std=gnu++11 -Wno-delete-non-virtual-dtor -Wall -Wextra")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wstrict-aliasing -pedantic -fmax-errors=5 -Werror=return-type")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wunreachable-code -Wcast-align -Wcast-qual")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel")
        set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option")
    else()
        set(CMAKE_CXX_FLAGS                "-std=gnu++11 -Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type")
    endif()

    set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0 -DDEBUG ")
    set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -DNDEBUG")


    if(WARNINGS_AS_ERRORS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werr")
    endif()

    if(OPT_ASAN)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
    else()
        if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING AND NOT WIN32 AND CODE_COVERAGE)
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
        endif()
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

        if(NOT EXISTS ${CMAKE_CXX_COMPILER})
            message( FATAL_ERROR "Clang++ not found. " )
        endif()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register") # -Qunused-arguments
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-const-variable -Wno-unused-parameter")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    endif()
endif()

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()
