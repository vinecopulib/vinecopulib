if (WIN32)
    set(gtest_disable_pthreads on) #needed in MinGW
    set(gtest_force_shared_crt ON CACHE BOOL "Always use msvcrt.dll" FORCE)
endif ()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include_directories(SYSTEM ${external_includes})
include_directories(${CMAKE_BINARY_DIR}/generated/test/)
add_subdirectory(src_test)

foreach (app ${unit_tests})
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} gtest vinecopulib ${CMAKE_THREAD_LIBS_INIT} ${test_libs} src_test)
endforeach ()
