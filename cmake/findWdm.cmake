include(cmake/DownloadProject.cmake)
download_project(
        PROJ wdm
        GIT_REPOSITORY https://github.com/tnagler/wdm.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET 1
)
set(wdm_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/wdm-src/include")
file(GLOB_RECURSE wdm_main ${wdm_INCLUDE_DIRS}/wdm.hpp)
file(GLOB_RECURSE wdm_hpp ${wdm_INCLUDE_DIRS}/wdm/*.hpp)
