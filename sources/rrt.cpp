#include "rrt.h"
#include "geometry.h"
#include <algorithm>
#include <iostream>


// find pointer to nearest point in tree to sample

tree_node <double>* find_nearest(
    std::vector <tree_node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance
) {
    tree_node <double> *nearest(&tree[0]);

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
    const tree_node <double>& A, 
    const point <double>& B, 
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
) {

    point <double> new_cfg(new_state(A, B, stepsize, collision_check));
    return new_cfg == B;
}


// find nearest tree_node we can connect directly to some point

tree_node <double>* connect_closest (
    std::vector <tree_node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
) {
    tree_node <double> *nearest(nullptr);

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

// get new state, along line from nearest tree_node to sample tree_node

point <double> new_state (
    const tree_node <double>& nearest, 
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


/*
Function which implements RRT algorithm.

start, goal are the terminal configurations on the path
joint_limits contains the min/max values of the joint variables
joint_wraps tells us if the joint variable wraps around (has circular topology)
collision_check returns if configuration is in free space
get_sample returns the next sample
num_samples is the number of samples used in algorithm
*/


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
    if (collision_check(start)) {
        throw std::domain_error("Start has to be in free space");
    }
    if (collision_check(goal)) {
        throw std::domain_error("Goal has to be in free space");
    }

    // initialize tree with starting point
    std::vector <tree_node <double>> tree;
    tree.reserve(num_samples + 2);
    tree.push_back(tree_node <double>(start, nullptr));

    // std::cout << "collision_check(" << start << ") = " << collision_check(start) << std::endl;
    // return output_rrt(false, tree);

    for (size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
        // get new (random) sample
        point <double> sample(get_sample());

        // find nearest point in RRT
        tree_node <double> *nearest(find_nearest(tree, sample, distance));

        // find new configuration
        point <double> new_cfg;
        try {
            new_cfg = new_state(*nearest, sample, stepsize, collision_check);
        }
        catch (std::length_error err) {
            throw;
        }

        // add new configuration to tree
        tree.push_back(tree_node <double>(new_cfg, nearest));
    }

    // find nearest tree_node we can connect to directly from goal
    tree_node <double>* nearest(connect_closest(tree, goal, distance, stepsize, collision_check));

    if (nearest == nullptr) {
        // output failed, return tree with false flag
        return output_rrt(false, tree);
    }

    // add this connection to tree
    tree.push_back(tree_node <double>(goal, nearest));

    // return successful output
    return output_rrt(true, tree);
}