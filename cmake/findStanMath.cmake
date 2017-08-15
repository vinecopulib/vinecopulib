set(STAN_MATH "math-${STAN_MATH_VERSION}")

if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/math/stan/math.hpp
        OR NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes/include/cvodes/cvodes.h)
    set(STAN_MATH_GITHUB https://github.com/stan-dev/math)
    file(DOWNLOAD
            ${STAN_MATH_GITHUB}/archive/v${STAN_MATH_VERSION}.tar.gz
            ${PROJECT_BINARY_DIR}/${STAN_MATH}.tar.gz)

    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_BINARY_DIR}/${STAN_MATH}.tar.gz
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    file(COPY ${PROJECT_BINARY_DIR}/${STAN_MATH}/stan/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/math/stan)

    file(GLOB CVODES ${PROJECT_BINARY_DIR}/${STAN_MATH}/lib/*/include/cvodes/cvodes.h)
    string(REPLACE "/include/cvodes/cvodes.h" "" CVODES ${CVODES})
    file(COPY ${CVODES}/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes)

    if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/math/stan/math.hpp
            OR NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes/include/cvodes/cvodes.h)
        message(SEND_ERROR "Not able to find or download stan_math.")
    endif()
endif()

set(STAN_MATH_INCLUDE_DIRS
        ${PROJECT_BINARY_DIR}/generated/vinecopulib/math
        ${PROJECT_BINARY_DIR}/generated/vinecopulib/cvodes/include)