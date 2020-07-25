#include "fmt.h"
#include "geometry.h"
#include "planning.h"
#include <iostream>
#include <queue>


// implements connect procedure on fmt graph

std::vector <labeled_node <double>> induced_graph (
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
                vertex.add_neighbor(&other_vertex);
            }
        }
    }

    return lgraph;
}

// comparisson function for priority queue

bool compare::operator() (labeled_node <double> * const lhs, labeled_node <double> * const rhs) {
    return lhs->get_distance() > rhs->get_distance();
}


// reconstruct path from shortest path tree

std::vector <point <double>> reconstruct_path (
    labeled_node <double>* start, 
    labeled_node <double>* goal
) {
    std::vector <point<double>> path;
    labeled_node <double> *cur_node(goal);

    // iterate over path in reverse order
    while (cur_node) {
        path.push_back(cur_node->get_point());

        // adjust current node
        cur_node = cur_node->get_backpointer();
    }

    std::reverse(path.begin(), path.end());
    return path;
}


/*
Function which implements FMT algorithm.

start, goal are the terminal configurations on the path
joint_limits contains the min/max values of the joint variables
joint_wraps tells us if the joint variable wraps around (has circular topology)
collision_check returns if configuration is in free space
get_sample returns the next sample
num_samples is the number of samples used in algorithm
*/

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
) {
    // obtain all samples
    std::vector <point <double>> samples(get_samples(num_samples, get_sample, collision_check));

    size_t free_samples = samples.size();
    std::vector <node <double>> graph(free_samples + 2);

    // add start node to graph
    graph[0] = node <double>(start);

    // add free samples to graph
    for (size_t i = 0; i < free_samples; ++i) {
        graph[i + 1] = node <double>(samples[i]);
    }

    // add goal node to graph
    graph[free_samples + 1] = node <double>(goal);

    // connect nodes amongst themselves
    std::vector <labeled_node <double>> lgraph = induced_graph(graph, distance, collision_check, radius, stepsize);
    std::vector <std::pair <point<double>, point<double>>> elist;

    std::cout << std::endl << "Induced graph connected " << "(" << lgraph.size() << " nodes)" << std::endl;

    // unmarked nodes are treated as infinitely distant
    // add start node to V_open
    std::priority_queue <labeled_node <double>*, std::vector <labeled_node<double>*>, compare> pq;
    pq.push(&lgraph[0]);
    lgraph[0].add_mark();
    bool found_path = false;

    while (!pq.empty()) {
        labeled_node <double>* z = pq.top();
        pq.pop();

        // if we´ve reached goal
        if (z == &lgraph.back()) {
            found_path = true;
            break;
        }

        // for all neighbors of z
        for (auto&& x : z->get_unmarked_neighbors()) {
            labeled_node <double>* y_near = nullptr;

            // get marked neighbors of x in open set, find the closest one
            for (auto&& y : x->get_marked_neighbors()) {
                if (y_near == nullptr || 
                    y->get_distance() + distance(y->get_point(), x->get_point()) 
                    < y_near->get_distance() + distance(y_near->get_point(), x->get_point())) {

                    y_near = y;
                }
            }

            // if path from y_near to x is clear, we´ve found the shortest path to x
            if (y_near != nullptr && path_clear(y_near->get_point(), x->get_point(), collision_check, stepsize)) {
                x->set_distance(y_near->get_distance() + distance(y_near->get_point(), x->get_point()));
                x->add_mark();
                x->set_backpointer(y_near);
                pq.push(x);
                elist.push_back(std::make_pair(x->get_point(), y_near->get_point()));
            }
        }

        // add z to V_closed
        z->close();
        z->remove_mark();
    }

    std::vector <std::vector <point <double>>> paths;

    // if there is no path return empty path
    if (!found_path) {
        paths.resize(1);
        return output(paths, elist);
    }

    // reconstruct path in treeified graph
    std::vector <point <double>> path = reconstruct_path(&lgraph[0], &lgraph.back());
    paths.push_back(path);
    return output(paths, elist);
}