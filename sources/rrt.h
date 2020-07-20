#pragma once

#include <vector>
#include <functional>
#include "geometry.h"
#include "planning.h"


// Function prototypes

tree_node <double>* find_nearest (
    std::vector <tree_node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance
);

tree_node <double>* connect_closest (
    std::vector <tree_node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check
);

bool is_reachable(
    const tree_node <double>& A, 
    const point <double>& B, 
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
);

tree_node <double>* connect_closest (
    std::vector <tree_node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
);

point <double> new_state (
    const tree_node <double>& nearest, 
    const point <double>& sample, 
    double stepsize, 
    const std::function <bool(point<double>)>& collision_check,
    const double epsilon = 1e-6
);

std::vector <std::pair<point<double>, point<double>>> edgelist (std::vector <tree_node <double>>& tree);

output rrt (
    const point <double> start,
    const point <double> goal,
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <bool(point<double>)>& collision_check,
    const std::function <point<double>()>& get_sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const unsigned int num_samples,
    const double stepsize
);