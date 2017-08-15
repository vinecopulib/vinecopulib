set(STAN "stan-${STAN_VERSION}")

if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan/stan/version.hpp)
    set(STAN_GITHUB https://github.com/stan-dev/stan)
    file(DOWNLOAD
            ${STAN_GITHUB}/archive/v${STAN_VERSION}.tar.gz
            ${PROJECT_BINARY_DIR}/${STAN}.tar.gz)

    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_BINARY_DIR}/${STAN}.tar.gz
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    file(COPY ${PROJECT_BINARY_DIR}/${STAN}/src/stan/
            DESTINATION ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan/stan)

    if(NOT EXISTS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan/stan/version.hpp)
        message(SEND_ERROR "Not able to find or download stan.")
    endif()
endif()

set(STAN_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/generated/vinecopulib/stan)