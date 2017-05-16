find_package(NLopt REQUIRED)

# this is needed for GCC 6 because we later use these as SYSTEM headers
# see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=70129
#     https://issues.apache.org/jira/browse/THRIFT-3828
#     https://gitlab.kitware.com/cmake/cmake/issues/16291
get_filename_component(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIRS} REALPATH)
