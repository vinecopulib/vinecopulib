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

#include "include/structselect_tools.hpp"

namespace structselect_tools {
    
    //! Create base tree of the vine
    //! 
    //!  The base tree is a star on d + 1 variables, where the conditioning
    //!  set of each edge consists of a single number. When building the next
    //!  tree, the edges become vertices. Because the base graph was a star
    //!  all edges are allowed by the proximity condition, and the edges will
    //!  have a conditioning set consisting of the two vertex indices. This
    //!  will be the first actual tree of the vine.
    //!  
    //!  @param data nxd matrix of copula data.
    //!  @return A VineTree object containing the base graph.
    VineTree make_base_tree(const MatXd& data) {
        int d = data.cols();
        VineTree base_tree(d);
        // a star connects the root node (d) with all other nodes
        for (int target = 0; target < d; ++target) {
            // add edge and extract edge iterator
            auto e = add_edge(d, target, base_tree).first;
            
            // add edge data & info
            base_tree[e].hfunc1 = data.col(boost::target(e, base_tree));

            std::vector<int> conditioning(1);
            conditioning[0] = boost::target(e, base_tree);
            base_tree[e].conditioning = conditioning;

            std::vector<int> conditioned(0);
            conditioned.reserve(d - 2);
            base_tree[e].conditioned = conditioned;
        }

        return base_tree;
    }
    
    //! Build next tree of the vine
    //! 
    //! The next tree is found the following way:
    //!     1. Edges of the previous tree become edges in the new tree.
    //!     2. All edges allowed by the proximity condition are added to the new
    //!        graph.
    //!     3. Collapse the new graph to a maximum spanning tree for edge weight
    //!        |tau|.
    //! 
    //! @param prev_tree tree T_{k}.
    //! @param tree T_{k+1}.
    VineTree build_next_tree(VineTree& prev_tree) 
    {
        // edges of old graph become vertices in new graph
        auto new_tree = edges_as_vertices(prev_tree);
                
        return new_tree;
    }
    
    //! Convert edge set into vertex set of a new graph
    //! 
    //! Further information about the structure is passed along:
    //!     - conditioning/conditioned set,
    //!     - indices of vertices connected by the edge in the previous tree.
    //! 
    //! @param tree T_{k}.
    //! @param A edge-less graph of vertices, each representing one edge of the
    //! previous tree.
    VineTree edges_as_vertices(const VineTree& prev_tree) {
        // start with full graph
        int d = num_edges(prev_tree);
        VineTree new_tree(d);
        
        // cut & paste information from previous tree
        int i = 0;
        for (auto e : boost::edges(prev_tree)) {
            new_tree[i].hfunc1 = prev_tree[e].hfunc1;
            new_tree[i].hfunc2 = prev_tree[e].hfunc2;
            new_tree[i].conditioning = prev_tree[e].conditioning;
            new_tree[i].conditioned = prev_tree[e].conditioned;
            new_tree[i].prev_edge_indices.reserve(2);
            new_tree[i].prev_edge_indices.push_back(boost::source(e, prev_tree));
            new_tree[i].prev_edge_indices.push_back(boost::target(e, prev_tree));
            ++i;
        }
        
        return new_tree;
    }
    
}
