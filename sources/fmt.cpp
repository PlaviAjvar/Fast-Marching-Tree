#include "fmt.h"
#include "geometry.h"
#include <algorithm>
#include <iostream>

/*
Function which implements FMT algorithm.

Function returns vector of joint space points, along the path. 

start, goal are the terminal configurations on the path
joint_limits contains the min/max values of the joint variables
joint_wraps tells us if the joint variable wraps around (has circular topology)
collision_check returns if configuration is in free space
get_sample returns the next sample
num_samples is the number of samples used in algorithm
*/


// returns path from start to goal in tree
// the endpoints are the first and last entry of vector

std::vector <point<double>> output::get_path(
    std::vector <node<double>>& tree
) {
    std::vector <point<double>> path;
    node <double> *cur_node(&tree.back());

    // iterate over path in reverse order
    while (cur_node != nullptr) {
        path.push_back(cur_node->get_point());

        // adjust current node
        cur_node = cur_node->get_parent();
    }

    std::reverse(path.begin(), path.end());
    return path;
}



// get edgelist from node-base graph

std::vector <std::pair<point<double>, point<double>>> output::edgelist_from_tree(
    std::vector <node <double>>& tree
) {
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


// find pointer to nearest point in tree to sample

node <double>* find_nearest(
    std::vector <node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance
) {
    node <double> *nearest(&tree[0]);

    for (auto&& vertex : tree) {
        if (distance(vertex.get_point(), sample) < distance(nearest->get_point(), sample)) {
            nearest = &vertex;
        }
    }

    // std::cout << nearest << " " << nearest->get_point() << std::endl;
    return nearest;
}


// check if point is reachable from another point

bool is_reachable(
    const node <double>& A, 
    const point <double>& B, 
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
) {

    point <double> new_cfg(new_state(A, B, stepsize, collision_check));
    return new_cfg == B;
}


// find nearest node we can connect directly to some point

node <double>* connect_closest (
    std::vector <node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
) {
    node <double> *nearest(nullptr);

    for (auto&& vertex : tree) {
        if (is_reachable(vertex, sample, stepsize, collision_check)) {
            if (nearest == nullptr || distance(vertex.get_point(), sample) < distance(nearest->get_point(), sample)) {
                nearest = &vertex;
            }
        }
    }

    // std::cout << nearest << " " << nearest->get_point() << std::endl;
    return nearest;
}

// walk along direction from nearest to sample

point <double> walk (
    const point <double>& nearest, 
    const point <double>& sample, 
    const double scaler
) {

    point <double> diff(sample - nearest);
    return nearest + diff * scaler;
}


// get new state, along line from nearest node to sample node

point <double> new_state (
    const node <double>& nearest, 
    const point <double>& sample, 
    double stepsize, 
    const std::function <bool(point<double>)>& collision_check,
    const double epsilon
) {

    point <double> new_cfg(nearest.get_point());

    for (double scaler = stepsize; scaler <= 1 + epsilon; scaler += stepsize) {
        // point 1 step further down line
        point <double> step_further = walk(nearest.get_point(), sample, scaler);

        // if no collision update new configuration
        if (collision_check(step_further)) {
            return new_cfg;
        }
        
        new_cfg = step_further;
    }
    
    return sample;
}

// Helper function for comparing solution quality
// implements RRT algorithm, same function template

output rrt(
    const point <double> start,
    const point <double> goal,
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <bool(point<double>)>& collision_check,
    const std::function <point<double>()>& get_sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const unsigned int num_samples,
    const double stepsize
) {
    // initialize tree with starting point
    std::vector <node <double>> tree;
    tree.reserve(num_samples + 2);
    tree.push_back(node <double>(start, nullptr));

    for (size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
        // get new (random) sample
        point <double> sample(get_sample());

        // find nearest point in RRT
        node <double> *nearest(find_nearest(tree, sample, distance));

        // find new configuration
        point <double> new_cfg;
        try {
            new_cfg = new_state(*nearest, sample, stepsize, collision_check);
        }
        catch (std::length_error err) {
            std::cout << err.what() << "\n";
            throw;
        }

        // add new configuration to tree
        tree.push_back(node <double>(new_cfg, nearest));
    }

    // find nearest node we can connect to directly from goal
    node <double>* nearest(connect_closest(tree, goal, distance, stepsize, collision_check));

    if (nearest == nullptr) {
        // output failed, return tree with false flag
        return output(false, tree);
    }

    // add this connection to tree
    tree.push_back(node <double>(goal, nearest));

    // return successful output
    return output(true, tree);
}