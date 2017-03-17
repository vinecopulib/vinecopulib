#pragma once

// differentiate between compilation of the library and of the client code
#if defined (_WIN32)
#    if defined(vinecopulib_EXPORTS)
#        define  VINECOPULIB_EXPORT __declspec(dllexport)
#    else
#        define  VINECOPULIB_EXPORT __declspec(dllimport)
#    endif
#else
#    define VINECOPULIB_EXPORT
#endif
