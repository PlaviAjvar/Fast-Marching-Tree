#include "matplotlibcpp.h"
#include "geometry.h"
#include "fmt.h"
#include <functional>
namespace plt = matplotlibcpp;

class test {

public:
    test () {}

    test (const test& P) {}

    bool test_collisionA(point <double> P, const double eps = 1e-3) const {
        bool col = (P[0] < eps || P[0] > 1-eps || P[1] < eps || P[1] > 1-eps) || (P[0] < 0.6 && (P[0] < 0.2 + P[1] && P[0] > -0.2 + P[1]));
        // std::cout << P << " : ";

        // if (col) std::cout << "True";
        // else std::cout << "False";
        
        // std::cout << std::endl;
        return col;
    }

    point <double> get_sample2d() const {
        double x = double(rand()) / RAND_MAX;
        double y = double(rand()) / RAND_MAX;
        point <double> sample(std::vector<double>{x, y});
        return sample;
    }
};


int main () {
    srand(time(NULL));

    // simple collision test
    double x_start, y_start;
    double x_goal, y_goal;
    x_start = 0.1;
    y_start = 0.4;
    x_goal = 0.4;
    y_goal = 0.1;

    point <double> start(std::vector<double>{x_start, y_start});
    point <double> goal(std::vector<double>{x_goal, y_goal});
    std::vector <std::pair<double,double>> joint_limits{{0,1}, {0,1}};
    unsigned int num_samples = 500;
    double step_size = 1e-2;

    test testA;
    std::function<bool(point<double>)> test_collision = [&testA] (point <double> P) { return testA.test_collisionA(P); };
    std::function<point<double>()> get_sample = [&testA] () { return testA.get_sample2d(); };
    std::function <double(const point <double>&, const point <double>&)> distance(euclidean_distance <double>);

    output result = rrt (
        start, goal, joint_limits, 
        test_collision, get_sample,
        distance, num_samples, step_size
    );

    std::cout << "In main" << std::endl;
    std::vector <double> xv, yv, xp, yp;

    for (const auto& vertex : result.get_edgelist()) {
        double x_par = (vertex.first)[0];
        double y_par = (vertex.first)[1];
        double x = (vertex.second)[0];
        double y = (vertex.second)[1];

        xv.push_back(x);
        yv.push_back(y);
        xp.push_back(x_par);
        yp.push_back(y_par);
    }

    std::cout << "Tree has " << xv.size() << " nodes" << std::endl;

    if (result.succeeded()) {
        double x_par, y_par;
        x_par = x_start;
        y_par = y_start;

        for (const auto& vertex : result.get_path()) {
            double x = vertex[0];
            double y = vertex[1];
            xv.push_back(x);
            yv.push_back(y);
            xp.push_back(x_par);
            yp.push_back(y_par);
            x_par = x;
            y_par = y;
        }
    }
    else {
        std::cout << "No path found." << std::endl;
    }

    plt::xlim(0, 1);
    plt::ylim(0, 1);

    for (size_t i = 0; i < xv.size(); ++i) {
        if (i < result.get_edgelist().size())
            plt::plot(std::vector <double>{xp[i], xv[i]}, std::vector <double>{yp[i], yv[i]}, "b-");
        else 
            plt::plot(std::vector <double>{xp[i], xv[i]}, std::vector <double>{yp[i], yv[i]}, "r-");
        plt::pause(1e-6);
    }

    plt::plot(xv, yv, "r-");
    plt::show();
    // plt::save("./output.png");
}