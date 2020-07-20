#include "matplotlibcpp.h"
#include "geometry.h"
#include "planning.h"
#include "rrt.h"
#include "prm.h"
#include "fmt.h"
#include <functional>
#include <random>
namespace plt = matplotlibcpp;

class test {
private:
    std::vector <point<double>> starts, goals;
    std::vector <std::pair <double,double>> joint_limits;
    unsigned int num_samples;
    double stepsize, radius;
    std::function<bool(point<double>)> test_collision;
    std::function<point<double>()> get_sample;
    std::function <double(const point <double>&, const point <double>&)> distance;

public:
    test () {}

    test (
        const std::vector <point<double>>& _starts,
        const std::vector <point<double>>& _goals,
        const std::vector <std::pair <double,double>>& _joint_limits,
        const std::function<bool(point<double>)> _test_collision,
        const std::function<point<double>()> _get_sample,
        std::function <double(const point <double>&, const point <double>&)> _distance,
        const unsigned int _num_samples,
        const double _stepsize,
        const double _radius
    ) : starts(_starts), goals(_goals), joint_limits(_joint_limits), num_samples(_num_samples), stepsize(_stepsize), radius(_radius),
    test_collision(_test_collision), get_sample(_get_sample), distance(_distance) {

        if (_starts.size() != _goals.size()) {
            throw std::length_error("Number of starting points has to equal number of goal points");
        }

        if (_joint_limits.size() != (_starts.begin())->get_dimension()) {
            throw std::domain_error("Number of joints has to equal dimension of configuration space points");
        }

        if (_stepsize < 0 || _radius < 0) {
            throw std::domain_error("Stepsize and radius have to be positive");
        }
    }

    test (
        const point <double>& start,
        const point <double>& goal,
        const std::vector <std::pair <double,double>>& _joint_limits,
        const std::function<bool(point<double>)> _test_collision,
        const std::function<point<double>()> _get_sample,
        std::function <double(const point <double>&, const point <double>&)> _distance,
        const unsigned int _num_samples,
        const double _stepsize,
        const double _radius
    ) : joint_limits(_joint_limits), num_samples(_num_samples), stepsize(_stepsize), radius(_radius),
    test_collision(_test_collision), get_sample(_get_sample), distance(_distance) {

        starts.push_back(start);
        goals.push_back(goal);

        if (_joint_limits.size() != (start.get_dimension())) {
            throw std::domain_error("Number of joints has to equal dimension of configuration space points");
        }

        if (_stepsize < 0 || _radius < 0) {
            throw std::domain_error("Stepsize and radius have to be positive");
        }
    }

    output run_test (const std::string& algorithm) const {
        if (algorithm != "PRM" && starts.size() > 1) {
            throw std::domain_error("Only PRM can have multiple start and goal points");
        }

        point <double> start = starts[0];
        point <double> goal = goals[0];
        output out;

        if (algorithm == "RRT") {
            out = rrt(start, goal, joint_limits, test_collision, get_sample, distance, num_samples, stepsize);
        }
        else if (algorithm == "FMT") {
            //out = fmt(start, goal, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else if (algorithm == "PRM") {
            out = prm(starts, goals, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else {
            throw std::domain_error("Algorithm label has to be either RRT, PRM or FMT");
        }

        return out;
    }
};


class test_battery {
public:
    // TEST A

    static point <double> startA () {
        return point <double>(std::vector<double>{0.1, 0.4});
    }

    static point <double> goalA () {
        return point <double>(std::vector<double>{0.4, 0.1});
    }

    static point <double> startA2 () {
        return point <double>(std::vector <double>{0.4, 0.8});
    }

    static point <double> goalA2 () {
        return point <double>(std::vector <double>{0.8, 0.9});
    }

    static std::vector <point <double>> startsA () {
        return std::vector <point <double>>{startA(), startA2()};
    }

    static std::vector <point<double>> goalsA () {
        return std::vector <point <double>>{goalA(), goalA2()};
    }

    static std::vector <std::pair<double,double>> joint_limitsA () {
        return std::vector <std::pair<double,double>>{{0,1}, {0,1}};
    }

    static std::function <double(const point <double>&, const point <double>&)> distanceA () {
        return euclidean_distance <double>;
    }

    constexpr static unsigned int num_samplesA = 500;
    constexpr static double stepsizeA = 1e-2;
    constexpr static double radiusA = 1e-2;

    static bool test_collisionA(point <double> P) {
        const double eps = 1e-3;
        bool col = (P[0] < eps || P[0] > 1-eps || P[1] < eps || P[1] > 1-eps) || (P[0] < 0.6 && (P[0] < 0.2 + P[1] && P[0] > -0.2 + P[1]));
        return col;
    }

    static size_t add_obstacle_edgesA (
        std::vector <point2d <double>>& current, 
        std::vector <point2d <double>>& previous, 
        const std::function<bool(point<double>)> test_collision
    ) {
        previous.push_back(point2d <double>(0, 0.2));
        current.push_back(point2d <double>(0.6, 0.8));

        previous.push_back(point2d <double>(0.6, 0.8));
        current.push_back(point2d <double>(0.6, 0.4));

        previous.push_back(point2d <double>(0.6, 0.4));
        current.push_back(point2d <double>(0.2, 0));

        return 3;
    }

    static point <double> get_sampleA () {
        return get_sample(joint_limitsA(), 2);
    }

    // TEST B

    static point <double> get_sample(const std::vector <std::pair <double,double>>& joint_limits, size_t dimension) {
        std::vector <double> components;
        std::random_device r;

        for (size_t dim = 0; dim < dimension; ++dim) {
            double lower = joint_limits[dim].first;
            double upper = joint_limits[dim].second;
            std::uniform_real_distribution<double> unif(lower, upper);
            std::default_random_engine re(r());
            components.push_back(unif(re));
        }

        point <double> sample(components);
        return sample;
    }
};

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
    const std::function <size_t(std::vector <point2d <double>>&, std::vector <point2d <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
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
            plt::pause(0.001);
        }
    }

    if (!save_image) {
        plt::show();
    }
    else {
        plt::save("./output.png");
    }
}


int main (int argc, char *argv[]) {
    test_battery tb;
    test testA_sq, testA_mq, testB_sq, testB_mq;

    try {
        testA_sq = test(tb.startA(), tb.goalA(), tb.joint_limitsA(), tb.test_collisionA, tb.get_sampleA, tb.distanceA(), tb.num_samplesA, tb.stepsizeA, tb.radiusA);
        testA_mq = test(tb.startsA(), tb.goalsA(), tb.joint_limitsA(), tb.test_collisionA, tb.get_sampleA, tb.distanceA(), tb.num_samplesA, tb.stepsizeA, tb.radiusA);
    }
    catch (std::logic_error err) {
        std::cout << err.what() << std::endl;
        throw;
    }

    std::string algorithm = "RRT";
    std::string test_label = "A";
    bool one_by_one = false;

    for (size_t i = 0; i < argc; ++i) {
        std::string flag(argv[i]);
        if (flag == "RRT" || flag == "PRM" || flag == "FMT") {
            algorithm = flag;
        }
        if (flag == "A" || flag == "B") {
            test_label = flag;
        }
        if (flag == "-obo") {
            one_by_one = true;
        }
    }

    output out;

    try {
        if (test_label == "A") {
            if (algorithm == "PRM") out = testA_mq.run_test(algorithm);
            else out = testA_sq.run_test(algorithm);
        }
        // else if (test_label == "B") {
        //     if (algorithm == "PRM") out = testB_mq.run_test(algorithm);
        //     else out = testB_sq.run_test(algorithm);
        // }
    }
    catch (std::logic_error err) {
        std::cout << err.what() << std::endl;
        throw;
    }

    // plot graph
    if (test_label == "A") plot_graph(out, tb.test_collisionA, tb.add_obstacle_edgesA, one_by_one);
    //else if (test_label == "B")) plot_graph(out, tb.test_collisionB, (algorithm != "PRM"));

    return 0;
}