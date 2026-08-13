# Compiling the documentation snippets is what keeps the examples in docs/*.dox
# from drifting: every fenced block on those pages is a \snippet of a file here.
file(GLOB snippet_sources CONFIGURE_DEPENDS *.cpp)

add_custom_target(doc_snippets)
foreach(snippet ${snippet_sources})
    get_filename_component(name ${snippet} NAME_WE)
    add_executable(doc_snippet_${name} EXCLUDE_FROM_ALL ${snippet})
    target_link_libraries(doc_snippet_${name} PRIVATE vinecopulib)
    add_dependencies(doc_snippets doc_snippet_${name})
endforeach()
