#include "prm.h"
#include "geometry.h"
#include <queue>
#include <iostream>

// implements connect procedure on prm graph

std::vector <labeled_node <double>> prm_connect(
    const std::vector <node <double>>& graph, 
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check,
    const double radius,
    const double stepsize
) {
    std::vector <labeled_node <double>> lgraph(graph.begin(), graph.end());

    for (auto&& vertex : lgraph) {
        for (auto&& other_vertex : lgraph) {
            // if they are different points and distance is smaller than radius
            if (vertex.get_point() != other_vertex.get_point() && distance(vertex.get_point(), other_vertex.get_point()) < radius) {
                // check if there is a collision
                if (path_clear(vertex.get_point(), other_vertex.get_point(), collision_check, stepsize)) {
                    vertex.add_neighbor(&other_vertex);
                }
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
    labeled_node <double>* goal
) {

    // BFS starting from start, ending at goal

    std::queue <labeled_node <double>*> q;
    q.push(start);
    start->add_mark();

    while (!q.empty()) {
        labeled_node <double>* front_node = q.front();
        q.pop();    

        // if goal is reached break bfs
        if (front_node == goal) {
            break;
        }

        for (auto&& neighbor : front_node->get_labeled_neighbors()) {
            if (!(neighbor->is_marked())) {
                neighbor->add_mark();
                neighbor->set_backpointer(front_node);
                q.push(neighbor);
            }
        }
    }

    if (!goal->is_marked()) return std::vector <point <double>>();

    return backtrack(start, goal);
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
        std::vector <point <double>> path (find_path(&lgraph[free_samples + 2 * query_idx], &lgraph[free_samples + 2 * query_idx + 1]));

        for (auto&& vertex : lgraph) {
            vertex.remove_mark();
        }

        paths.push_back(path);
    }

    return output(paths, lgraph);
}