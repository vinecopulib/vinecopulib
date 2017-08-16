set(STAN "stan-${STAN_VERSION}")
set(STAN_MATH "math-${STAN_MATH_VERSION}")

if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan/version.hpp)
    if(NOT EXISTS ${PROJECT_BINARY_DIR}/${STAN}.tar.gz)
        set(STAN_GITHUB https://github.com/stan-dev/stan)
        file(DOWNLOAD
                ${STAN_GITHUB}/archive/v${STAN_VERSION}.tar.gz
                ${PROJECT_BINARY_DIR}/${STAN}.tar.gz)
    endif()

    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_BINARY_DIR}/${STAN}.tar.gz
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    file(COPY ${PROJECT_BINARY_DIR}/${STAN}/src/stan/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan
            FILES_MATCHING PATTERN "*.hpp")
    file(COPY ${PROJECT_BINARY_DIR}/${STAN}/src/stan/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan
            FILES_MATCHING PATTERN "*.h")

    if (BUILD_TESTING)
        file(COPY ${PROJECT_BINARY_DIR}/${STAN}/src/test/test-models/good/optimization/rosenbrock.stan
                DESTINATION ${PROJECT_BINARY_DIR}/bin)
    endif()

    if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan/version.hpp)
        message(SEND_ERROR "Not able to find or download stan.")
    endif()
endif()

if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan/math.hpp)

    if(NOT EXISTS ${PROJECT_BINARY_DIR}/${STAN_MATH}.tar.gz)
        set(STAN_MATH_GITHUB https://github.com/stan-dev/math)
        file(DOWNLOAD
                ${STAN_MATH_GITHUB}/archive/v${STAN_MATH_VERSION}.tar.gz
                ${PROJECT_BINARY_DIR}/${STAN_MATH}.tar.gz)
    endif()

    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_BINARY_DIR}/${STAN_MATH}.tar.gz
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    file(COPY ${PROJECT_BINARY_DIR}/${STAN_MATH}/stan/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan)

    file(GLOB CVODES ${PROJECT_BINARY_DIR}/${STAN_MATH}/lib/*/include/cvodes/cvodes.h)
    string(REPLACE "include/cvodes/cvodes.h" "" CVODES ${CVODES})
    file(COPY ${CVODES}/include/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes_headers)
    #string(REPLACE "cvodes/cvodes.h" "" CVODES ${CVODES})
    #file(COPY ${CVODES}
    #        DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes_headers,
     #       FILES_MATCHING PATTERN "*.h")

    if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers/stan/math.hpp)
        message(SEND_ERROR "Not able to find or download stan_math.")
    endif()
endif()

file(GLOB_RECURSE stan_sources_inst
        ${PROJECT_BINARY_DIR}/${STAN}/src/stan/lang/*_inst.cpp)
file(GLOB_RECURSE stan_sources_def
        ${PROJECT_BINARY_DIR}/${STAN}/src/stan/lang/*_def.cpp)

#set(stanc_source ${PROJECT_BINARY_DIR}/${STAN}/src/test/test-models/stanc.cpp)
set(STAN_INCLUDE_DIRS
        ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan_headers
        ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes_headers)
