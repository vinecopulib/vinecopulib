# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  add_custom_target(git_revision.hpp ALL
    git log -1 "--format=format:#define GIT_REVISION \"%H\"%n" HEAD > include/git_revision.hpp
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM
  )
endif()

install(
  FILES ${CMAKE_SOURCE_DIR}/include/git_revision.hpp
  DESTINATION "${include_install_dir}"
)
