#include "matplotlibcpp.h"
#include "geometry.h"
#include "planning.h"
#include "rrt.h"
#include "prm.h"
#include "fmt.h"
#include <functional>
namespace plt = matplotlibcpp;

class test {

public:
    test () {}

    test (const test& P) {}

    bool test_collisionA(point <double> P, const double eps = 1e-3) const {
        bool col = (P[0] < eps || P[0] > 1-eps || P[1] < eps || P[1] > 1-eps) || (P[0] < 0.6 && (P[0] < 0.2 + P[1] && P[0] > -0.2 + P[1]));
        return col;
    }

    point <double> get_sample2d() const {
        double x = double(rand()) / RAND_MAX;
        double y = double(rand()) / RAND_MAX;
        point <double> sample(std::vector<double>{x, y});
        return sample;
    }
};


size_t add_obstacle_edges(
    std::vector <point2d <double>>& current, 
    std::vector <point2d <double>>& previous, 
    const std::function<bool(point<double>)>& test_collision
) {

    // dummy code, insert real code

    previous.push_back(point2d <double>(0, 0.2));
    current.push_back(point2d <double>(0.6, 0.8));

    previous.push_back(point2d <double>(0.6, 0.8));
    current.push_back(point2d <double>(0.6, 0.4));

    previous.push_back(point2d <double>(0.6, 0.4));
    current.push_back(point2d <double>(0.2, 0));

    return 3;
}


size_t add_graph_edges(
    std::vector <point2d <double>>& current, 
    std::vector <point2d <double>>& previous, 
    const std::vector <std::pair<point<double>, point<double>>>& elist
) {

    for (const auto& vertex : elist) {
        point2d <double> prev_point(vertex.first);
        point2d <double> cur_point(vertex.second);

        previous.push_back(prev_point);
        current.push_back(cur_point);
    }

    return current.size();
}


size_t add_path_edges(
    std::vector <point2d <double>>& current, 
    std::vector <point2d <double>>& previous, 
    const std::vector <std::vector <point <double>>>& paths
) {

    for (auto&& path : paths) {
        for (size_t i = 1; i < path.size(); ++i) {
            previous.push_back(path[i-1]);
            current.push_back(path[i]);
        }
    }

    return current.size();
}


void plot_graph(
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    bool delay_active = true,
    bool save_image = false
) {

    std::vector <point2d <double>> current, previous;

    // add obstacle edges
    size_t cutoffA = add_obstacle_edges(current, previous, test_collision);

    // add graph edges
    size_t cutoffB = add_graph_edges(current, previous, result.get_edgelist());

    // add path edges
    size_t cutoffC = add_path_edges(current, previous, result.get_paths());

    plt::xlim(0, 1);
    plt::ylim(0, 1);

    for (size_t i = 0; i < cutoffC; ++i) {
        std::string lineflag;

        if (i < cutoffA) lineflag = "g-";
        else if (i < cutoffB) lineflag = "b-";
        else lineflag = "r-";

        plt::plot(std::vector <double>{previous[i].getx(), current[i].getx()}, std::vector <double>{previous[i].gety(), current[i].gety()}, lineflag); 
        
        if (delay_active) {
            plt::pause(0.01);
        }
    }

    if (!save_image) {
        plt::show();
    }
    else {
        plt::save("./output.png");
    }
}


int main () {
    srand(time(NULL));

    // simple collision test
    double x_start, y_start;
    double x_goal, y_goal;
    x_start = 0.1;
    y_start = 0.4;
    x_goal = 0.4;
    y_goal = 0.1;

    point <double> start2(std::vector <double>{0.4, 0.8});
    point <double> goal2(std::vector <double>{0.8, 0.9});

    point <double> start(std::vector<double>{x_start, y_start});
    point <double> goal(std::vector<double>{x_goal, y_goal});
    std::vector <point <double>> starts{start, start2};
    std::vector <point <double>> goals{goal, goal2};
    std::vector <std::pair<double,double>> joint_limits{{0,1}, {0,1}};
    unsigned int num_samples = 500;
    double step_size = 1e-2;
    double radius = 1e-2;

    test testA;
    std::function<bool(point<double>)> test_collision = [&testA] (point <double> P) { return testA.test_collisionA(P); };
    std::function<point<double>()> get_sample = [&testA] () { return testA.get_sample2d(); };
    std::function <double(const point <double>&, const point <double>&)> distance(euclidean_distance <double>);

    output_rrt resultRRT = rrt (
        start, goal, joint_limits, 
        test_collision, get_sample,
        distance, num_samples, step_size
    );

    output resultPRM = prm (
        starts, goals, joint_limits,
        test_collision, get_sample,
        distance, num_samples, step_size, radius
    );

    std::cout << "In main" << std::endl;

    // plot graph
    // plot_graph(resultRRT, test_collision);
    plot_graph(resultPRM, test_collision, false);

    return 0;
}