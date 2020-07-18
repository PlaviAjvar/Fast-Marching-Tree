#pragma once

#include "geometry.h"
#include <vector>

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
public:
    tree_node () {}

    tree_node (const tree_node <real>& vertex) {
        node<real>::sample = vertex.get_point();
        node<real>::neighbor = vertex.get_neighbors();
    }

    tree_node (
        const point <real> _sample,
        tree_node <real> *parent
    ) {

        node<real>::sample = _sample;
        node<real>::neighbor.resize(1);
        node<real>::neighbor[0] = parent;
    }

    tree_node <real>* get_parent () const {
        if (node<real>::neighbor.size() == 0) return nullptr;
        return dynamic_cast <tree_node <real>*> (node<real>::neighbor[0]);
    }
};


template <typename real>
class labeled_node : public node <real> {
private:
    bool mark;
    labeled_node <real>* backpointer;

public:
    labeled_node () {}

    labeled_node (const node <real>& vertex) {
        node<real>::sample = vertex.get_point();
        node<real>::neighbor = vertex.get_neighbors();
    }

    labeled_node (const point <real> _sample) {
        node<real>::sample = _sample;
    }

    labeled_node (
        const point <real> _sample, 
        std::vector <node<real>*> _neighbor
    ) : node<real>(_sample, _neighbor) {}

    std::vector <labeled_node<real>*> get_labeled_neighbors () const {
        std::vector <labeled_node <real>*> labeled_neighbor;
        labeled_neighbor.reserve(node<real>::neighbor.size());

        for (auto&& vertex : node<real>::neighbor) {
            labeled_neighbor.push_back(dynamic_cast <labeled_node <real>*>(vertex));
        }

        return labeled_neighbor;
    }

    labeled_node <double>* get_backpointer() const {
        return backpointer;
    }

    bool is_marked () const {
        return mark;
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
    output () {}

    output (const output& other) : path(other.path), edge_list(other.edge_list) {}

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