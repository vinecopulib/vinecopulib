/*
Copyright 2016 Thibault Vatter, Thomas Nagler

This file is part of vinecopulib.

vinecopulib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vinecopulib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vinecopulib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VINECOPULIB_STRUCTSELECT_TOOLS_HPP
#define VINECOPULIB_STRUCTSELECT_TOOLS_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/graph_utility.hpp>

// to allow for (auto e : boost::edges(g)) notation
namespace std
{
    template <class T>
    T begin(const std::pair<T,T>& eItPair) { return eItPair.first; }
    
    template <class T>
    T end(const std::pair<T,T>& eItPair) { return eItPair.second; }
}

namespace structselect_tools {
    
}


#endif
