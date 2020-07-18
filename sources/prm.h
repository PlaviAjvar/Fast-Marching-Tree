#pragma once

#include "geometry.h"
#include "planning.h"

output prm (
    const std::vector <point <double>>& start,
    const std::vector <point <double>>& goal,
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <bool(point<double>)>& collision_check,
    const std::function <point<double>()>& get_sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const unsigned int num_samples,
    const double stepsize,
    const double radius
);

std::vector <labeled_node <double>> prm_connect(
    const std::vector <node <double>>& graph, 
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check,
    const double radius,
    const double stepsize
);

std::vector <point <double>> find_path (
    labeled_node <double>* start,
    labeled_node <double>* goal
);

std::vector <point<double>> backtrack (
    labeled_node <double>* start,
    labeled_node <double>* goal
);

bool path_clear(
    const point <double> A,
    const point <double> B,
    const std::function <bool(point<double>)>& collision_check,
    const double stepsize,
    const double epsilon = 1e-6
);