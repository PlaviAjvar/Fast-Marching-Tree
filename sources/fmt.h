#pragma once

#include <vector>
#include <functional>
#include "geometry.h"

// defines node class for graphs

template <typename real>
class node {
private:
    point <real> sample;
    node <real> *parent;

public:
    node () {}

    node (const point <real> sample, node <real> * const parent) {
        this->sample = sample;
        this->parent = parent;
    }

    node (const node<real>& vertex) : sample(vertex.get_point()), parent(vertex.get_parent()) {}

    point <real> get_point () const {
        return sample;
    }

    node <real>* get_parent () const {
        return parent;
    }
};


// class representing output from algorithms
class output {
private:
    bool success;
    std::vector <point <double>> path;
    std::vector <std::pair<point<double>, point<double>>> edge_list;

    std::vector <std::pair<point<double>, point<double>>> edgelist_from_tree(
        std::vector <node <double>>& tree
    );

    std::vector <point<double>> get_path(
        std::vector <node <double>>& tree
    );

public:
    output () {}

    output (
        bool _success, 
        std::vector <point<double>> _path, 
        std::vector <std::pair<point<double>, point<double>>> elist
    ) : success(_success), path(_path), edge_list(elist) {}

    output (const output& other) : success(other.success), path(other.path), edge_list(other.edge_list) {}

    // if we have vector of nodes as input
    output (
        bool _success, 
        std::vector <node <double>>& tree
    ) : success(_success) {
        
        path = get_path(tree);
        edge_list = edgelist_from_tree(tree);
    }

    bool succeeded () const {
        return success;
    }

    std::vector <point <double>> get_path () const {
        return path;
    }

    std::vector <std::pair<point<double>, point<double>>> get_edgelist () const {
        return edge_list;
    }
};


// Function prototypes

std::vector <point<double>> get_path (std::vector <node<double>>& tree);

node <double>* find_nearest (
    std::vector <node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance
);

node <double>* connect_closest (
    std::vector <node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check
);

bool is_reachable(
    const node <double>& A, 
    const point <double>& B, 
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
);

node <double>* connect_closest (
    std::vector <node <double>>& tree,
    const point <double>& sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const double stepsize,
    const std::function <bool(point<double>)>& collision_check
);

point <double> walk (
    const point <double>& nearest, 
    const point <double>& sample, 
    const double scaler
);

point <double> new_state (
    const node <double>& nearest, 
    const point <double>& sample, 
    double stepsize, 
    const std::function <bool(point<double>)>& collision_check,
    const double epsilon = 1e-6
);

std::vector <std::pair<point<double>, point<double>>> edgelist (std::vector <node <double>>& tree);

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