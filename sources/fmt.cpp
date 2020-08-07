#include "fmt.h"
#include "geometry.h"
#include "planning.h"
#include <iostream>
#include <queue>
#include <algorithm>
#include <iterator>


// implements connect procedure on fmt graph

std::vector <labeled_node <double>> induced_graph (
    const std::vector <node <double>>& graph, 
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const std::function <bool(point<double>)>& collision_check,
    const double radius,
    const double stepsize
) {
    std::vector <labeled_node <double>> lgraph(graph.begin(), graph.end());

    for (auto it = lgraph.begin(); it != lgraph.end(); ++it) {
        for (auto jt = std::next(it); jt != lgraph.end(); ++jt) {
            if (distance(it->get_point(), jt->get_point()) < radius) {
                it->add_neighbor(&(*jt));
                jt->add_neighbor(&(*it));
            }
        }
    }

    return lgraph;
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
    if (collision_check(start)) {
        throw std::domain_error("Start has to be in free space");
    }
    if (collision_check(goal)) {
        throw std::domain_error("Goal has to be in free space");
    }
    
    // obtain all samples
    std::vector <point <double>> samples(get_samples(num_samples, get_sample, collision_check));

    std::cout << "num_samples = " << samples.size() << std::endl;
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

    std::cout << "Induced graph connected " << "(" << lgraph.size() << " nodes)" << std::endl << std::endl;

    // unmarked nodes are treated as infinitely distant
    // add start node to V_open
    std::priority_queue <labeled_node <double>*, std::vector <labeled_node<double>*>, compare<double>> pq;
    pq.push(&lgraph[0]);
    lgraph[0].add_mark();
    lgraph[0].set_distance(0);
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


// function for obtaining radius for FMT* algorithm
// as explained in the original paper
// gamma is a tuning parameter, which has to satisfy relation given in paper

double fmtradius (
    const unsigned int num_samples,
    const unsigned int dimension,
    const double gamma
) {
    double nf = num_samples;
    double d = dimension;
    double radius = gamma * pow(log(nf) / nf, 1 / d);
    return radius;
}


// Lebesgue meassure of configuration space
// since configuration space is N-d box
// mu = dim(1) * dim(2) * ... * dim(d)

double lebesgue (
    const std::vector <std::pair<double,double>>& joint_limits
) {
    double mu = 1;
    for (const auto& lim : joint_limits) {
        mu *= std::abs(lim.second - lim.first);
    }
    return mu;
}

// volume of d-dimensional unit hypersphere
// volume = pi^(n/2) * R^n / gamma(n/2 + 1)

double volume (
    const unsigned int dim
) {
    double d = dim;
    return pow(M_PI, d / 2) / tgamma(1 + d / 2);
}


// asume whole configuration space is free space (worst case case)

double fmtgamma (
    const std::vector <std::pair<double,double>>& joint_limits,
    const double eta
) {
    double d = joint_limits.size();
    double gamma = 2 * (1 + eta) * pow(1 / d, 1 / d) * pow(lebesgue(joint_limits) / volume(d), 1 / d);
    return gamma;
}


// algorithm implementing generic FMT*, with self-tuning radius parameter
output fmtstar (
    const point <double> start,
    const point <double> goal,
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <bool(point<double>)>& collision_check,
    const std::function <point<double>()>& get_sample,
    const std::function <double(const point <double>&, const point <double>&)> distance,
    const unsigned int num_samples,
    const double stepsize,
    const double eta
) {
    if (eta < 0) {
        std::cout << "Warning: Eta has to be greater than zero for asympototic optimality guarantee." << std::endl << std::endl;
    }

    double gamma = fmtgamma(joint_limits, eta);
    double radius = fmtradius(num_samples, joint_limits.size(), gamma);
    double adjust = adjust_weight(joint_limits, distance);  // adjust for weighed euclidean metric (if it´s active)
    radius *= adjust;

    std::cout << std::endl << "Radius(" << num_samples << ") = " << radius << std::endl;
    std::cout << "Adjust(weighed euclidean) = " << adjust << std::endl << std::endl;
    return fmt(start, goal, joint_limits, collision_check, get_sample, distance, num_samples, stepsize, radius);
}

double adjust_weight(
    const std::vector <std::pair<double,double>>& joint_limits,
    const std::function <double(const point <double>&, const point <double>&)> distance
) {

    double scale = 1;
    std::vector <double> components(joint_limits.size());

    for (size_t dim = 0; dim < joint_limits.size(); ++dim) {
        double dif = joint_limits[dim].second - joint_limits[dim].first;

        // set vector with leftmost corner point
        for (size_t dim = 0; dim < joint_limits.size(); ++dim) {
            components[dim] = joint_limits[dim].first;
        }
        point <double> lower(components);

        // change current coordinate to max limit
        components[dim] = joint_limits[dim].second;
        point <double> upper(components);

        // calculate distance and compare to coordinate difference
        double dist = distance(lower, upper);
        double ratio = dist / dif;

        // adjust for current dimension scaling (multiply by sqrt(ratio))
        scale *= sqrt(ratio);
    }

    // apply pow(scale, 1/d) to obtain actual contribution
    double d = joint_limits.size();
    return pow(scale, 1 / d);
}