// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

namespace vinecopulib {

namespace tools_os {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        static const std::string rm="del ";
        static const std::string stanc="stanc.exe ";
#else
        static const std::string rm="rm ";
        static const std::string stanc="stanc ";
#endif
}

}
