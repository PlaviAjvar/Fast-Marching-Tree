#pragma once

#include "geometry.h"
#include <vector>
#include <iostream>

// defines node class for graphs

template <typename real>
class node {
protected:
    point <real> sample;
    std::vector <node <real>*> neighbor;

public:
    node () {}

    node (const point <real> _sample) : sample(_sample) {}

    node (
        const point <real> _sample, 
        std::vector <node<real>*> _neighbor
    ) : sample(_sample), neighbor(_neighbor) {}

    node (const node<real>& vertex) : sample(vertex.get_point()), neighbor(vertex.get_neighbors()) {}

    point <real> get_point () const {
        return sample;
    }

    std::vector <node<real>*> get_neighbors () const {
        return neighbor;
    }

    void add_neighbor (node <real> * const vertex) {
        neighbor.push_back(vertex);
    }

    virtual ~node() = default;
};


template <typename real>
class tree_node : public node <real> {
private:
    real distance;

public:
    tree_node () {}

    tree_node (const tree_node <real>& vertex) {
        node<real>::sample = vertex.get_point();
        node<real>::neighbor = vertex.get_neighbors();
        distance = vertex.get_distance();
    }

    tree_node (
        const point <real> _sample,
        tree_node <real> *parent,
        const real dist = 0
    ) {

        node<real>::sample = _sample;
        node<real>::neighbor.resize(1);
        node<real>::neighbor[0] = parent;
        distance = dist;
    }

    tree_node <real>* get_parent () const {
        if (node<real>::neighbor.size() == 0) return nullptr;
        return dynamic_cast <tree_node <real>*> (node<real>::neighbor[0]);
    }

    real get_distance () const {
        return distance;
    }
};


template <typename real>
class labeled_node : public node <real> {
private:
    bool mark;
    bool closed;
    double distance;
    labeled_node <real>* backpointer;
    const double inf = 1e10;

public:
    labeled_node () {
        backpointer = nullptr;
        mark = false;
        closed = false;
        distance = inf;
    }

    labeled_node (const node <real>& vertex) {
        node<real>::sample = vertex.get_point();
        node<real>::neighbor = vertex.get_neighbors();
        backpointer = nullptr;
        mark = false;
        closed = false;
        distance = inf;
    }

    labeled_node (const point <real> _sample) {
        node<real>::sample = _sample;
        backpointer = nullptr;
        mark = false;
        closed = false;
        distance = inf;
    }

    labeled_node (
        const point <real> _sample, 
        std::vector <node<real>*> _neighbor
    ) : node<real>(_sample, _neighbor) {
        
        backpointer = nullptr;
        mark = false;
        closed = false;
        distance = inf;
    }

    std::vector <labeled_node<real>*> get_labeled_neighbors () const {
        std::vector <labeled_node <real>*> labeled_neighbors;
        labeled_neighbors.reserve(node<real>::neighbor.size());

        for (auto&& vertex : node<real>::neighbor) {
            labeled_neighbors.push_back(dynamic_cast <labeled_node <real>*>(vertex));
        }

        return labeled_neighbors;
    }

    std::vector <labeled_node<real>*> get_marked_neighbors () const {
        std::vector <labeled_node <real>*> marked_neighbors;
        marked_neighbors.reserve(node <real>::neighbor.size());

        for (auto&& vertex : node<real>::neighbor) {
            labeled_node <real>* lvertex = dynamic_cast <labeled_node <real>*>(vertex);
            if (lvertex -> is_marked()) {
                marked_neighbors.push_back(lvertex);
            }
        }

        return marked_neighbors;
    }

    std::vector <labeled_node<real>*> get_unmarked_neighbors () const {
        std::vector <labeled_node <real>*> unmarked_neighbors;
        unmarked_neighbors.reserve(node <real>::neighbor.size());

        for (auto&& vertex : node<real>::neighbor) {
            labeled_node <real>* lvertex = dynamic_cast <labeled_node <real>*>(vertex);
            if (!lvertex->is_marked() && !lvertex->is_closed()) {
                unmarked_neighbors.push_back(lvertex);
            }
        }

        return unmarked_neighbors;
    }
    
    labeled_node <double>* get_backpointer() const {
        return backpointer;
    }

    bool is_marked () const {
        return mark;
    }

    bool is_closed () const {
        return closed;
    }

    bool is_inf (const double eps = 1) const {
        return std::abs(distance - inf) < eps;
    }

    double get_distance () const {
        return distance;
    }

    void set_backpointer (labeled_node <double> * const bp) {
        backpointer = bp;
    }

    void add_mark () {
        mark = true;
    }

    void remove_mark () {
        mark = false;
    }

    void close() {
        closed = true;
    }

    void set_distance (const double dist) {
        distance = dist;
    }

    void clear () {
        mark = false;
        closed = false;
        distance = inf;
        backpointer = nullptr;
    }
};


// class representing output from algorithms
class output {
protected:
    std::vector <std::vector <point <double>>> path;
    std::vector <std::pair<point<double>, point<double>>> edge_list;

    std::vector <std::pair<point<double>, point<double>>> edgelist_from_graph (
        std::vector <labeled_node <double>>& graph
    ) const;

public:
    output () {
        path.resize(1);
    }

    output (const output& other) : path(other.get_paths()), edge_list(other.get_edgelist()) {}

    output (
        const std::vector <std::vector <point<double>>>& _path, 
        const std::vector <std::pair<point<double>, point<double>>>& elist
    ) : path(_path), edge_list(elist) {}

    output (
        const std::vector <std::vector <point<double>>>& _path, 
        std::vector <labeled_node <double>>& graph
    ) : path(_path) {

        edge_list = edgelist_from_graph(graph);
    }

    std::vector <std::vector <point <double>>> get_paths () const {
        return path;
    }

    std::vector <std::pair<point<double>, point<double>>> get_edgelist () const {
        return edge_list;
    }

    void add_path (const std::vector <point <double>>& new_path) {
        path.push_back(new_path);
    }
};


// output for the RRT algorithm

class output_rrt : public output {
private:
    std::vector <std::pair<point<double>, point<double>>> edgelist_from_tree(
        std::vector <tree_node <double>>& tree
    ) const;

    std::vector <point<double>> get_path(
        std::vector <tree_node <double>>& tree
    ) const;

public:
    output_rrt () {}

    output_rrt (const output_rrt& other) : output(other.get_paths(), get_edgelist()) {}

    output_rrt (
        std::vector <point<double>> _path, 
        std::vector <std::pair<point<double>, point<double>>> elist
    ) {

        edge_list = elist;
        path.resize(1);
        path[0] = _path;
    }

    // if we have vector of nodes as input
    output_rrt (
        bool _success, 
        std::vector <tree_node <double>>& tree
    ) {
        if (_success) {
            path.resize(1);
            path[0] = get_path(tree);
        }
        edge_list = edgelist_from_tree(tree);
    }

    bool succeeded () const {
        return !path[0].empty();
    }
};

std::vector <point <double>> get_samples (
    unsigned int num_samples,
    const std::function <point<double>()>& get_sample,
    const std::function <bool(point<double>)>& collision_check,
    const size_t max_iter = 100000000
);

template <typename real>
class compare {
public:
// comparisson function for priority queue

    bool operator() (labeled_node <real> * const lhs, labeled_node <real> * const rhs) {
        return lhs->get_distance() > rhs->get_distance();
    }
};

