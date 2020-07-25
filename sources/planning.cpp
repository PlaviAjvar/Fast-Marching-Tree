#include "planning.h"
#include <iostream>

// returns path from start to goal in tree
// the endpoints are the first and last entry of vector

std::vector <point<double>> output_rrt::get_path(
    std::vector <tree_node<double>>& tree
) const {
    std::vector <point<double>> path;
    tree_node <double> *cur_node(&tree.back());

    // iterate over path in reverse order
    while (cur_node != nullptr) {
        path.push_back(cur_node->get_point());

        // adjust current node
        cur_node = cur_node->get_parent();
    }

    std::reverse(path.begin(), path.end());
    return path;
}



// get edgelist from node-based graph

std::vector <std::pair<point<double>, point<double>>> output_rrt::edgelist_from_tree (
    std::vector <tree_node <double>>& tree
) const {
    std::vector <std::pair<point<double>, point<double>>> elist;
    elist.reserve(tree.size());

    for (auto&& node : tree) {
        // std::cout << "child = " << node.get_point() << std::endl;
        point <double> parent, child;
        child = node.get_point();
        if (node.get_parent() != nullptr) {
            parent = (node.get_parent())->get_point();
        }
        else {
            parent = node.get_point();
        }
        // std::cout << "parent = " << parent << std::endl;
        elist.push_back({parent, child});
    }

    return elist;    
}

template <class T>
std::ostream& operator<< (std::ostream& os, std::vector <T> v) {
    for(auto&& e : v) std::cout << e << " ";
    return os;
}

// obtain edgelist from graph

std::vector <std::pair<point<double>, point<double>>> output::edgelist_from_graph (
        std::vector <labeled_node <double>>& graph
) const {
    std::vector <std::pair<point<double>, point<double>>> edge_list;

    for (auto&& vertex : graph) {
        for (auto&& neighbor : vertex.get_neighbors()) {
            edge_list.push_back(std::make_pair(vertex.get_point(), neighbor->get_point()));
        }
    }

    return edge_list;
}


// helper function for obtaining list of samples

std::vector <point <double>> get_samples(
    unsigned int num_samples,
    const std::function <point<double>()>& get_sample,
    const std::function <bool(point<double>)>& collision_check
) {

    // obtain all samples
    std::vector <point <double>> samples;
    samples.reserve(num_samples);

    for (size_t i = 0; i < num_samples; ++i) {
        point <double> sample = get_sample();
        
        // if sample is in free space add it
        if (!collision_check(sample)) {
            samples.push_back(sample);
        }
    }

    return samples;
}