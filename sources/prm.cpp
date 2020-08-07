#include "prm.h"
#include "geometry.h"
#include <queue>
#include <iostream>
#include <algorithm>

// implements connect procedure on prm graph

std::vector <labeled_node <double>> prm_connect(
    const std::vector <node <double>>& graph, 
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check,
    const double radius,
    const double stepsize
) {
    std::vector <labeled_node <double>> lgraph(graph.begin(), graph.end());

    for (auto it = lgraph.begin(); it != lgraph.end(); ++it) {
        for (auto jt = std::next(it); jt != lgraph.end(); ++jt) {
            // check if there is a collision
            if (distance(it->get_point(), jt->get_point()) < radius && path_clear(it->get_point(), jt->get_point(), collision_check, stepsize)) {
                it->add_neighbor(&(*jt));
                jt->add_neighbor(&(*it));
            }
        }
    }

    return lgraph;
}


// backtrack solution after bfs

std::vector <point<double>> backtrack (
    labeled_node <double>* start,
    labeled_node <double>* goal
) {

    labeled_node <double>* cur_node = goal;
    std::vector <point <double>> path;

    while (cur_node != start) {
        path.push_back(cur_node->get_point());
        cur_node = cur_node->get_backpointer();
    }

    path.push_back(start->get_point());
    std::reverse(path.begin(), path.end());
    return path;
}


// find path from start_point to goal_point through graph

std::vector <point <double>> find_path (
    labeled_node <double>* start,
    labeled_node <double>* goal,
    const std::function <double(const point <double>&, const point <double>&)> distance
) {

    // Apply Dijkstra´s algorithm to find shortest path through graph

    std::priority_queue <labeled_node <double>*, std::vector <labeled_node <double>*>, compare<double>> pq;
    pq.push(start);
    start->set_distance(0);

    while (!pq.empty()) {
        labeled_node <double>* top_node = pq.top();
        pq.pop();    

        // if node has already been processed
        if (top_node->is_marked()) {
            continue;
        }
        
        top_node->add_mark();

        // if goal is reached return solution
        if (top_node == goal) {
            return backtrack(start, goal);
        }

        // there is no path in this case
        if (top_node->is_inf()) {
            break;
        }

        for (auto&& neighbor : top_node->get_unmarked_neighbors()) {
            if (top_node->get_distance() + distance(top_node->get_point(), neighbor->get_point()) < neighbor->get_distance()) {
                neighbor->set_distance(top_node->get_distance() + distance(top_node->get_point(), neighbor->get_point()));
                neighbor->set_backpointer(top_node);
                pq.push(neighbor);
            }
        }
    }

    // if we´ve reached here a path has not been found
    return std::vector <point <double>>();
}

/*
Function which implements PRM algorithm.

start, goal are the terminal configurations on the path
joint_limits contains the min/max values of the joint variables
joint_wraps tells us if the joint variable wraps around (has circular topology)
collision_check returns if configuration is in free space
get_sample returns the next sample
num_samples is the number of samples used in algorithm
*/

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
) { 
    for (auto &st : start) {
        if (collision_check(st)) {
            throw std::domain_error("Starts have to be in free space");
        }
    }
    for (auto &go : goal) {
        if (collision_check(go)) {
            throw std::domain_error("Goals have to be in free space");
        }
    }

    std::vector <point <double>> samples(get_samples(num_samples, get_sample, collision_check));

    size_t free_samples = samples.size();
    size_t num_queries = start.size();
    std::vector <node <double>> graph(free_samples + 2 * num_queries);

    // add free samples to graph
    for (size_t i = 0; i < free_samples; ++i) {
        graph[i] = node <double>(samples[i]);
    }

    // add start and end nodes to graph
    for (size_t i = 0; i < num_queries; ++i) {
        graph[free_samples + 2*i] = node <double>(start[i]);
        graph[free_samples + 2*i + 1] = node <double>(goal[i]);
    }

    // connect nodes amongst themselves
    std::vector <labeled_node <double>> lgraph = prm_connect(graph, distance, collision_check, radius, stepsize);
    std::vector <std::vector <point <double>>> paths;

    std::cout << "PRM connected" << std::endl;

    // process queries
    for (size_t query_idx = 0; query_idx < num_queries; ++query_idx) {
        std::vector <point <double>> path (find_path(&lgraph[free_samples + 2 * query_idx], &lgraph[free_samples + 2 * query_idx + 1], distance));
        paths.push_back(path);
        
        for (auto&& vertex : lgraph) {
            vertex.clear();
        }
    }

    return output(paths, lgraph);
}