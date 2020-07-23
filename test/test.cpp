#include "matplotlibcpp.h"
#include "geometry.h"
#include "planning.h"
#include "rrt.h"
#include "prm.h"
#include "fmt.h"
#include <functional>
#include <random>
#include <sstream>
#include <cmath>
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
            out = fmt(start, goal, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else if (algorithm == "PRM") {
            out = prm(starts, goals, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else {
            throw std::domain_error("Algorithm label has to be either RRT, PRM or FMT");
        }

        // std::cout << start << " " << goal << std::endl;
       
        // for (auto&& edge : out.get_edgelist()) std::cout << edge.first << " -> " << edge.second << std::endl;

        return out;
    }
};


class test_battery {
public:
    /***************************
    TEST A (pointlike robot 2D)
    ***************************/

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
    constexpr static double radiusA = 1e-2; // 5 / num_samples

    static bool test_collisionA (const point <double> P) {
        const double eps = 1e-3;
        bool col = (P[0] < eps || P[0] > 1-eps || P[1] < eps || P[1] > 1-eps) || (P[0] < 0.6 && (P[0] < 0.2 + P[1] && P[0] > -0.2 + P[1]));
        return col;
    }

    static size_t add_obstacle_edgesA (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous, 
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

    /***************************
    TEST B (pointlike robot 2D)
    ***************************/

    static point <double> startB () {
        return point <double>(std::vector<double>{0.1, 0.3});
    }

    static point <double> goalB () {
        return point <double>(std::vector<double>{0.8, 0.8});
    }

    static point <double> startB2 () {
        return point <double>(std::vector <double>{0.1, 0.8});
    }

    static point <double> goalB2 () {
        return point <double>(std::vector <double>{0.8, 0.1});
    }

    static std::vector <point <double>> startsB () {
        return std::vector <point <double>>{startB(), startB2()};
    }

    static std::vector <point<double>> goalsB () {
        return std::vector <point <double>>{goalB(), goalB2()};
    }

    static std::vector <std::pair<double,double>> joint_limitsB () {
        return std::vector <std::pair<double,double>>{{0,1}, {0,1}};
    }

    static std::function <double(const point <double>&, const point <double>&)> distanceB () {
        return euclidean_distance <double>;
    }

    constexpr static unsigned int num_samplesB = 500;
    constexpr static double stepsizeB = 1e-2;
    constexpr static double radiusB = 1e-2; // 5 / num_samples


    static bool test_collisionB (const point <double> P) {
        const double eps = 1e-3;

        if (P[0] < eps || P[0] > 1-eps || P[1] < eps || P[1] > 1-eps) {
            return true;
        }
        if (P[0] < 0.2 && P[1] > 0.45 && P[1] < 0.55) {
            return true;
        }

        if (P[1] < 0.5 && P[0] < 0.6 && P[1] < P[0]) {
            return true;
        }

        if (P[0] > 0.7 && P[1] > 0.3 && P[1] < 0.7) {
            return true;
        }

        if (P[0] < 0.6 && P[0] > 0.5 && P[1] > 0.7) {
            return true;
        }

        return false;
    }

    static size_t add_obstacle_edgesB (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous, 
        const std::function<bool(point<double>)> test_collision
    ) {
        // left part
        previous.push_back(point2d <double>(0, 0.45));
        current.push_back(point2d <double>(0.2, 0.45));

        previous.push_back(point2d <double>(0.2, 0.45));
        current.push_back(point2d <double>(0.2, 0.55));

        previous.push_back(point2d <double>(0.2, 0.55));
        current.push_back(point2d <double>(0, 0.55));

        // bottom part
        previous.push_back(point2d <double>(0, 0));
        current.push_back(point2d <double>(0.5, 0.5));

        previous.push_back(point2d <double>(0.5, 0.5));
        current.push_back(point2d <double>(0.6, 0.5));

        previous.push_back(point2d <double>(0.6, 0.5));
        current.push_back(point2d <double>(0.6, 0));

        // right part
        previous.push_back(point2d <double>(1, 0.3));
        current.push_back(point2d <double>(0.7, 0.3));

        previous.push_back(point2d <double>(0.7, 0.3));
        current.push_back(point2d <double>(0.7, 0.7));

        previous.push_back(point2d <double>(0.7, 0.7));
        current.push_back(point2d <double>(1, 0.7));

        // upper part
        previous.push_back(point2d <double>(0.5, 1));
        current.push_back(point2d <double>(0.5, 0.7));

        previous.push_back(point2d <double>(0.5, 0.7));
        current.push_back(point2d <double>(0.6, 0.7));

        previous.push_back(point2d <double>(0.6, 0.7));
        current.push_back(point2d <double>(0.6, 1));

        return 12;
    }

    static point <double> get_sampleB () {
        return get_sample(joint_limitsB(), 2);
    }

    /***************************
    TEST C (pointlike robot 3D)
    ***************************/

    static point <double> startC () {
        return point <double>(std::vector<double>{0.1, 0.5, 0.5});
    }

    static point <double> goalC () {
        return point <double>(std::vector<double>{0.5, 0.1, 0.9});
    }

    static std::vector<point <double>> startsC () {
        std::vector <point <double>> starts{startC()};
        return starts;
    }

    static std::vector <point <double>> goalsC () {
        std::vector <point <double>> goals{goalC()};
        return goals;
    }

    static std::vector <std::pair<double,double>> joint_limitsC () {
        return std::vector <std::pair<double,double>>{{0,1}, {0,1}, {0,1}};
    }

    static std::function <double(const point <double>&, const point <double>&)> distanceC () {
        return euclidean_distance <double>;
    }

    constexpr static unsigned int num_samplesC = 2500;
    constexpr static double stepsizeC = 1e-2;
    constexpr static double radiusC = 1e-2; // same as old

    static bool test_collisionC (const point <double> P) {
        const double eps = 1e-3;

        if (P[0] < eps || P[0] > 1-eps || P[1] < eps || P[1] > 1-eps || P[2] < eps || P[2] > 1-eps) {
            return true;
        }

        if (P[0] < 0.4 && P[1] > 0.4 && P[1] < 0.6) {
            return false;
        }

        if (P[0] > 0.6 && P[1] > 0.4 && P[1] < 0.6) {
            return false;
        }

        if (P[0] > 0.4 && P[0] < 0.6) {
            return false;
        }

        return true;
    }

    static point <double> get_sampleC () {
        return get_sample(joint_limitsC(), 3);
    }

    static size_t add_obstacle_edgesC (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous, 
        const std::function<bool(point<double>)> test_collision
    ) {
        // left side
        previous.push_back(point2d <double>(0, 0.4));
        current.push_back(point2d <double>(0.4, 0.4));

        previous.push_back(point2d <double>(0.4, 0.4));
        current.push_back(point2d <double>(0.4, 0));

        previous.push_back(point2d <double>(0, 0.6));
        current.push_back(point2d <double>(0.4, 0.6));
        
        previous.push_back(point2d <double>(0.4, 0.6));
        current.push_back(point2d <double>(0.4, 1));

        // right side

        previous.push_back(point2d <double>(1, 0.4));
        current.push_back(point2d <double>(0.6, 0.4));

        previous.push_back(point2d <double>(0.6, 0.4));
        current.push_back(point2d <double>(0.6, 0));

        previous.push_back(point2d <double>(1, 0.6));
        current.push_back(point2d <double>(0.6, 0.6));
        
        previous.push_back(point2d <double>(0.6, 0.6));
        current.push_back(point2d <double>(0.6, 1));

        return 8;
    }

    // function which generates random sample

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
    std::vector <point <double>>& current, 
    std::vector <point <double>>& previous, 
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
    std::vector <point <double>>& current, 
    std::vector <point <double>>& previous, 
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

std::ostream& operator<< (std::ostream& os, std::vector <std::vector <double>> v) {
    for (auto row : v) {
        for (auto el : v) {
            os << el << " ";
        }
        os << std::endl;
    }
    return os;
}


// plot 2D graph
void plot_graph (
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
    bool delay_active = true,
    bool save_image = false
) {

    std::vector <point <double>> current, previous;

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
        else if (i < cutoffB) lineflag = "-b";
        else lineflag = "r-";

        point2d <double> prev(previous[i]);
        point2d <double> cur(current[i]);
        plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, lineflag); 
        
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

// get as many as possible discrete levels of green-blue
std::vector <std::string> get_levels () {
    std::vector <std::string> levels;

    for (unsigned int lev = 1; lev < 256; ++lev) {
        std::stringstream ss;
        ss << std::hex << lev;
        std::string res = ss.str();

        for(auto &ch : res) {
            ch = toupper(ch);
        }
        if (res.size() == 1) res = "0" + res;

        levels.push_back("#00" + res + "FF");
    }

    return levels;
}

// get discrete level of blue from z coordinate
std::string scale_color (
    const std::vector <std::string>& levels,
    const double z,
    const double zlim
) {
    size_t levelcnt = levels.size();
    double quant = zlim / levelcnt;
    size_t lev = z / quant;
    return levels[lev];
}


// plot 3D graph
void plot_3d (
    const output& result,
    const std::function<bool(point<double>)>& test_collision,
    const std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges,
    bool delay_active = true,
    bool save_image = false,
    double zlim = 1
) {

    std::vector <point <double>> current, previous;

    // add obstacle edges
    size_t cutoffA = add_obstacle_edges(current, previous, test_collision);

    // add graph edges
    size_t cutoffB = add_graph_edges(current, previous, result.get_edgelist());

    // add path edges
    size_t cutoffC = add_path_edges(current, previous, result.get_paths());

    plt::xlim(0, 1);
    plt::ylim(0, 1);

    std::vector <std::string> levels(get_levels());

    for (size_t i = 0; i < cutoffC; ++i) {
        std::string lineflag;
        point3d <double> prev(previous[i]);
        point3d <double> cur(current[i]);

        if (i < cutoffA) {
            lineflag = "g-";
        }
        else if (i < cutoffB) {
            // if blue scale it
            lineflag = scale_color(levels, cur.getz(), zlim);
        }
        else {
            lineflag = "r-";
        }

        plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, lineflag); 
        
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
    test testA_sq, testA_mq, testB_sq, testB_mq, testC, testC_mq;

    try {
        testA_sq = test(tb.startA(), tb.goalA(), tb.joint_limitsA(), tb.test_collisionA, tb.get_sampleA, tb.distanceA(), tb.num_samplesA, tb.stepsizeA, tb.radiusA);
        testA_mq = test(tb.startsA(), tb.goalsA(), tb.joint_limitsA(), tb.test_collisionA, tb.get_sampleA, tb.distanceA(), tb.num_samplesA, tb.stepsizeA, tb.radiusA);
        testB_sq = test(tb.startB(), tb.goalB(), tb.joint_limitsB(), tb.test_collisionB, tb.get_sampleB, tb.distanceB(), tb.num_samplesB, tb.stepsizeB, tb.radiusB);
        testB_mq = test(tb.startsB(), tb.goalsB(), tb.joint_limitsB(), tb.test_collisionB, tb.get_sampleB, tb.distanceB(), tb.num_samplesB, tb.stepsizeB, tb.radiusB);

        testC = test(tb.startC(), tb.goalC(), tb.joint_limitsC(), tb.test_collisionC, tb.get_sampleC, tb.distanceC(), tb.num_samplesC, tb.stepsizeC, tb.radiusC);
        testC_mq = test(tb.startsC(), tb.goalsC(), tb.joint_limitsC(), tb.test_collisionC, tb.get_sampleC, tb.distanceC(), tb.num_samplesC, tb.stepsizeC, tb.radiusC);
    }
    catch (std::logic_error err) {
        std::cout << err.what() << std::endl;
        throw;
    }
    

    std::string algorithm = "FMT";
    std::string test_label = "C";
    bool one_by_one = false;
    bool write_to_file = false;

    for (size_t i = 0; i < argc; ++i) {
        std::string flag(argv[i]);
        if (flag == "RRT" || flag == "PRM" || flag == "FMT") {
            algorithm = flag;
        }
        if (flag == "A" || flag == "B" || flag == "C") {
            test_label = flag;
        }
        if (flag == "point2D-A") flag = "A";
        if (flag == "point2D-B") flag = "B";
        if (flag == "point3D") flag = "C";
        if (flag == "-seq") {
            one_by_one = true;
        }
        if (flag == "-file") {
            write_to_file = true;
        }
    }

    output out;

    try {
        if (test_label == "A") {
            if (algorithm == "PRM") out = testA_mq.run_test(algorithm);
            else out = testA_sq.run_test(algorithm);
        }
        else if (test_label == "B") {
            if (algorithm == "PRM") out = testB_mq.run_test(algorithm);
            else out = testB_sq.run_test(algorithm);
        }
        else if (test_label == "C") {
            if (algorithm == "PRM") out = testC_mq.run_test(algorithm);
            else out = testC.run_test(algorithm);
        }
    }
    catch (std::logic_error err) {
        std::cout << err.what() << std::endl;
        throw;
    }

    for (size_t qry = 0; qry < out.get_paths().size(); ++qry) {
        std::vector <point <double>> path = out.get_paths()[qry];
        if (path.size() == 0) {
            std::cout << "Path for query " << qry + 1 << " not found." << std::endl;
        }
        else {
            std::cout << "Path for query " << qry + 1 << " found." << std::endl;
        }
    }
    
    std::cout << out.get_paths().size() << std::endl;

    // plot graph
    if (test_label == "A") plot_graph(out, tb.test_collisionA, tb.add_obstacle_edgesA, one_by_one, write_to_file);
    else if (test_label == "B") plot_graph(out, tb.test_collisionB, tb.add_obstacle_edgesB, one_by_one, write_to_file);
    else if (test_label == "C") plot_3d(out, tb.test_collisionC, tb.add_obstacle_edgesC, one_by_one, write_to_file, 1);

    return 0;
}