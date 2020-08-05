#pragma once

#include "geometry.h"
#include "planning.h"

output fmt (
    const point <double> start,
    const point <double> goal,
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <bool(point<double>)>& collision_check,
    const std::function <point<double>()>& get_sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const unsigned int num_samples,
    const double stepsize,
    const double radius
);

output fmtstar (
    const point <double> start,
    const point <double> goal,
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <bool(point<double>)>& collision_check,
    const std::function <point<double>()>& get_sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const unsigned int num_samples,
    const double stepsize,
    const double eta = 1
);

std::vector <labeled_node <double>> induced_graph (
    const std::vector <node <double>>& graph, 
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check,
    const double radius,
    const double stepsize
);

std::vector <point <double>> reconstruct_path (
    labeled_node <double>* start, 
    labeled_node <double*> goal
);

double fmtradius (
    const unsigned int num_samples,
    const unsigned int dimension,
    const double gamma
);

double lebesgue (
    const std::vector <std::pair<double,double>>& joint_limits
);

double volume (
    const unsigned int dim
);

double fmtgamma (
    const std::vector <std::pair<double,double>>& joint_limits,
    const double eta
);

double adjust_weight(
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <double(const point <double>&, const point <double>&)> distance
);