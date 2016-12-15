include_directories(${PROJECT_SOURCE_DIR}/)
include_directories(${PROJECT_SOURCE_DIR}/common)
#include_directories("../common")

include_directories(${external_includes})
#include_directories( ${OpenCV_INCLUDE_DIRS} )

# Create test executables
set(app_programs
        test_indep_bicop
        test_normal_bicop
        test_gumbel_bicop
        test_clayton_bicop)

foreach(app ${app_programs})
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app}
                          common
                          ${external_libs})
endforeach()
