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
#include <chrono>
#include <fstream>
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

        return out;
    }

    // some useful getters
    unsigned int get_num_samples () const {
        return num_samples;
    }

    double get_radius () const {
        return radius;
    }

    double get_stepsize() const {
        return stepsize;
    }
};


class test_battery {
public:
    // parameters
    constexpr static unsigned int num_samplesA = 500;
    constexpr static double stepsizeA = 1e-2;
    constexpr static double radiusA = 0.1; // 50 / num_samples

    constexpr static unsigned int num_samplesB = 500;
    constexpr static double stepsizeB = 1e-2;
    constexpr static double radiusB = 0.1; // 50 / num_samples

    constexpr static unsigned int num_samplesC = 1000;
    constexpr static double stepsizeC = 1e-2;
    constexpr static double radiusC = 1e-1; // same as old

    constexpr static unsigned int num_samplesD = 300; // since there are more collisions
    constexpr static double stepsizeD = 3e-2; // around PI * old_value
    constexpr static double radiusD = 3; // since weighed euclidean, have to scale
    // new_radius = total_length * strech_scaler * old_radius = 5 * 3 * 1e-1 = 3

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
        return get_sample_util(joint_limitsA(), 2);
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
        return get_sample_util(joint_limitsB(), 2);
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
        return get_sample_util(joint_limitsC(), 3);
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

    /***************************
    TEST D (2D planar manipulator)
    ***************************/

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesD_A () {
        // define points of obstacle
        point2d <double> A(5.1, 5);
        point2d <double> B(6.1, 4);
        point2d <double> C(9, 4);
        point2d <double> D(9, 7);
        point2d <double> E(5.1, 5);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesD_B () {
        // define points of obstacle
        point2d <double> A(-2, -2);
        point2d <double> B(0, -2);
        point2d <double> C(4, -6);
        point2d <double> D(-2, -8);
        point2d <double> E(-2, -2);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesD_C () {
        // define points of obstacle
        double h = 0.3;
        point2d <double> A(-2, 0);
        point2d <double> B(-h, 0);
        point2d <double> C(-h, 7);
        point2d <double> D(-2, 7);
        point2d <double> E(-2, 0);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesDT () {
        // define points of obstacle
        double h = 1;
        point2d <double> A(h, h);
        point2d <double> B(10, h);
        point2d <double> C(10, 10);
        point2d <double> D(h, 10);
        point2d <double> E(h, h);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static point2d <double> baseD () {
        point2d <double> base(0, 0);
        return base;
    }

    static std::vector <double> link_lengthsD () {
        std::vector <double> link_lengths{5, 5};
        return link_lengths;
    }

    static std::vector <std::pair <double,double>> joint_limitsD () {
        std::vector <std::pair <double,double>> lims {{-M_PI, M_PI}, {-M_PI, M_PI}};
        return lims;
    }

    static workspace <double> workspaceD (const std::string label) {
        // form polygons from points
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesA(edgesD_A());
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesB(edgesD_B());
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesC(edgesD_C());
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesT(edgesDT());
        std::vector <polygon <double>> obstacles;

        if (label == "D") {
            obstacles.push_back(polygon <double>(edgesA));
            obstacles.push_back(polygon <double>(edgesB));
            obstacles.push_back(polygon <double>(edgesC));
        }

        if (label == "DT") {
            obstacles.push_back(polygon <double>(edgesT));
        }

        // create 2D planar arm object
        point2d <double> base(baseD());
        std::vector <double> link_lengths{link_lengthsD()};
        std::vector <std::pair <double,double>> lims(joint_limitsD());

        // std::cout << "base = " << base << std::endl;
        // std::cout << "linklen = " << std::endl;
        // for (auto e : link_lengths) std::cout << e << ",";
        // std::cout << std::endl << "joint_lim = ";
        // for (auto e : lims) std::cout << "(" << e.first << " " << e.second << ")" << ",";
        // std::cout << std::endl;

        arm2d <double> planar(base, link_lengths, lims);

        workspace <double> ws(obstacles, planar);
        return ws;
    }

    static point <double> startD () {
        return point <double>(std::vector<double>{0, 0});
    }

    static point <double> goalD () {
        return point <double>(std::vector<double>{M_PI_2, -M_PI_2});
    }

    static point <double> goalDT () {
        return point <double>(std::vector <double>{-M_PI_2, -M_PI_2});
    }

    static std::vector<point <double>> startsD () {
        std::vector <point <double>> starts{startD()};
        return starts;
    }

    static std::vector <point <double>> goalsD () {
        std::vector <point <double>> goals{goalD()};
        return goals;
    }

    static std::vector <point <double>> goalsDT () {
        std::vector <point <double>> goals{goalDT()};
        return goals;
    }

    static double distanceD (const point <double>& A, const point <double>& B) {
        std::vector <double> linklen = workspaceD("D").get_link_lengths();
        return weighed_euclidean <double>(A, B, linklen);
    }

    static double distanceDT (const point <double>& A, const point <double>& B) {
        std::vector <double> linklen = workspaceD("DT").get_link_lengths();
        return weighed_euclidean <double>(A, B, linklen);
    }

    static bool test_collisionD (const point <double> config) {
        return workspaceD("D").collides(config);
    }

    static bool test_collisionDT (const point <double> config) {
        return workspaceD("DT").collides(config);
    }

    static point <double> get_sampleD () {
        return get_sample_util(joint_limitsD(), 2);
    }

    static size_t add_obstacle_edgesD (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous, 
        const std::function<bool(point<double>)> test_collision
    ) {
        size_t first = auto_add_edges(current, previous, edgesD_A());
        size_t second = auto_add_edges(current, previous, edgesD_B());
        return first + second;
    }

    // functions returning test info with test label as argument
    static point <double> start (std::string label) {
        if (label == "A") return startA();
        if (label == "B") return startB();
        if (label == "C") return startC();
        if (label == "D" || label == "DT") return startD();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static point <double> goal (std::string label) {
        if (label == "A") return goalA();
        if (label == "B") return goalB();
        if (label == "C") return goalC();
        if (label == "D") return goalD();
        if (label == "DT") return goalDT();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::vector<point <double>> starts (std::string label) {
        if (label == "A") return startsA();
        if (label == "B") return startsB();
        if (label == "C") return startsC();
        if (label == "D" || label == "DT") return startsD();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::vector <point <double>> goals (std::string label) {
        if (label == "A") return goalsA();
        if (label == "B") return goalsB();
        if (label == "C") return goalsC();
        if (label == "D") return goalsD();
        if (label == "DT") return goalsDT();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::vector <std::pair<double,double>> joint_limits (std::string label) {
        if (label == "A") return joint_limitsA();
        if (label == "B") return joint_limitsB();
        if (label == "C") return joint_limitsC();
        if (label == "D" || label == "DT") return joint_limitsD();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::function <double(const point <double>&, const point <double>&)> distance (std::string label) {
        if (label == "A" || label == "B" || label == "C") return euclidean_distance <double>;
        if (label == "D") return distanceD;
        if (label == "DT") return distanceDT;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static unsigned int num_samples (std::string label) {
        if (label == "A") return num_samplesA;
        if (label == "B") return num_samplesB;
        if (label == "C") return num_samplesC;
        if (label == "D" || label == "DT") return num_samplesD;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static double stepsize (std::string label) {
        if (label == "A") return stepsizeA;
        if (label == "B") return stepsizeB;
        if (label == "C") return stepsizeC;
        if (label == "D" || label == "DT") return stepsizeD;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static double radius (std::string label) {
        if (label == "A") return radiusA;
        if (label == "B") return radiusB;
        if (label == "C") return radiusC;
        if (label == "D" || label == "DT") return radiusD;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    std::function<bool(point<double>)> test_colision (std::string label) {
        if (label == "A") return test_collisionA;
        if (label == "B") return test_collisionB;
        if (label == "C") return test_collisionC;
        if (label == "D") return test_collisionD;
        if (label == "DT") return test_collisionDT;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    std::function<point<double>()> get_sample (std::string label) {
        if (label == "A") return get_sampleA;
        if (label == "B") return get_sampleB;
        if (label == "C") return get_sampleC;
        if (label == "D" || label == "DT") return get_sampleD;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    std::function <size_t(std::vector <point <double>>&, std::vector <point <double>>&, const std::function<bool(point<double>)>)> add_obstacle_edges (
        std::string label
    ) {
        if (label == "A") return add_obstacle_edgesA;
        if (label == "B") return add_obstacle_edgesB;
        if (label == "C") return add_obstacle_edgesC;
        if (label == "D" || label == "DT") return pass;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static workspace <double> getws (std::string label) {
        if (label == "D" || label == "DT") return workspaceD(label);
        
        // invalid test label
        throw std::domain_error("Invalid test label");
    }



    // get polygon edges from counterclockwise list of points
    static std::vector <std::pair <point2d<double>, point2d<double>>> edges_from_points (
        std::vector <point2d <double>> points
    ) {
        std::vector <std::pair <point2d<double>, point2d<double>>> edges(points.size() - 1);
        for (size_t i = 0; i < points.size()-1; ++i) {
            edges[i] = {points[i], points[i+1]};
        }
        return edges;
    }

    // function which generates random sample
    static point <double> get_sample_util (const std::vector <std::pair <double,double>>& joint_limits, size_t dimension) {
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

    // automatized adding of edges from edge list
    static size_t auto_add_edges (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous, 
        const std::vector <std::pair <point2d<double>, point2d<double>>>& edges
    ) {
        for (const auto& edge : edges) {
            previous.push_back(edge.first);
            current.push_back(edge.second);
        }
        return edges.size();
    }

    // function which passes
    static size_t pass (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous, 
        const std::function<bool(point<double>)> test_collision
    ) {
        return 0;
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
    const std::vector <std::pair <double, double>> joint_lims,
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

    plt::xlim(joint_lims[0].first, joint_lims[0].second);
    plt::ylim(joint_lims[1].first, joint_lims[1].second);

    for (size_t i = 0; i < cutoffC; ++i) {
        std::string lineflag;

        if (i < cutoffA) lineflag = "g-";
        else if (i < cutoffB) lineflag = "b-";
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
    const std::vector <std::pair <double, double>> joint_lims,
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

    plt::xlim(joint_lims[0].first, joint_lims[0].second);
    plt::ylim(joint_lims[1].first, joint_lims[1].second);

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

// display single snapshot
void disp_snapshot(
    const workspace <double>& ws, 
    const point <double> config, 
    const std::string color,
    const std::string name
) {
    auto robot = ws.get_robot();
    auto edges = robot.dir_kine(config);
    std::vector <double> xs;
    std::vector <double> ys;

    for (const auto& edge : edges) {
        xs.push_back(edge.first.getx());
        xs.push_back(edge.second.getx());
        ys.push_back(edge.first.gety());
        ys.push_back(edge.second.gety());
    }

    plt::named_plot(name, xs, ys, color + "-");
}


// display snapshots of 2D planar arm in its movement from start to goal
void display_snapshots (
    const workspace <double>& ws,
    const std::vector <point <double>>& path
) {
    std::vector <std::string> names{"start","2nd","3rd","4th","5th","goal"};
    std::vector <std::string> colors{"y", "m", "c", "r", "g", "b"};
    size_t num_snapshots = colors.size();

    for (const auto& polygon : ws.get_obstacles()) {
        for (const auto& edge : polygon.get_edges()) {
            point2d <double> prev = edge.first;
            point2d <double> cur = edge.second;
            plt::plot(std::vector <double>{prev.getx(), cur.getx()}, std::vector <double>{prev.gety(), cur.gety()}, "g-");
        }
    }

    disp_snapshot(ws, path[0], colors[0], names[0]);

    for (size_t i = 1; i < num_snapshots - 1; ++i) {
        double jump = double(path.size()) / num_snapshots;
        size_t path_idx = i * jump;
        disp_snapshot(ws, path[path_idx], colors[i], names[i]);
    }

    disp_snapshot(ws, path.back(), colors.back(), names.back());
    plt::legend();
    plt::show();
}


int main (int argc, char *argv[]) {
    test_battery tb;
    test test_sq, test_mq;

    std::string algorithm = "RRT";
    std::string test_label = "D";
    bool one_by_one = false;
    bool write_to_file = false;
    bool snapshot = false;
    std::string mode = "-normal";

    for (size_t i = 0; i < argc; ++i) {
        std::string flag(argv[i]);
        if (flag == "RRT" || flag == "PRM" || flag == "FMT") {
            algorithm = flag;
        }
        if (flag == "A" || flag == "B" || flag == "C" || flag == "D" || flag == "DT") {
            test_label = flag;
        }
        if (flag == "point2D-A") test_label = "A";
        if (flag == "point2D-B") test_label = "B";
        if (flag == "point3D") test_label = "C";
        if (flag == "arm2D" || flag == "planar") test_label = "D";
        if (flag == "trivial") test_label = "DT";
        if (flag == "-seq") {
            one_by_one = true;
        }
        if (flag == "-file") {
            write_to_file = true;
        }
        if (flag == "-sim" || flag == "-param" || flag == "-normal") {
            mode = flag;
        }
        if (flag == "-snapshot") {
            snapshot = true;
        }
    }

    // run algorithm normally and display results

    if (mode == "-normal") {
        auto begin = std::chrono::high_resolution_clock::now();
        output out;

        try {
            test_sq = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                tb.get_sample(test_label), tb.distance(test_label), tb.num_samples(test_label), tb.stepsize(test_label), tb.radius(test_label));

            out = test_sq.run_test(algorithm);
        }
        catch (std::logic_error err) {
            std::cout << err.what() << std::endl;
            throw;
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::cout << std::endl;

        for (size_t qry = 0; qry < out.get_paths().size(); ++qry) {
            std::vector <point <double>> path = out.get_paths()[qry];
            if (path.size() == 0) {
                std::cout << "Path for query " << qry + 1 << " not found" << std::endl;
            }
            else {
                std::cout << "Path for query " << qry + 1 << " found" << std::endl;
            }
        }

        std::cout << std::endl << "Elapsed time: " << elapsed.count() * 1e-9 << " seconds" << std::endl;

        double path_len;
        auto path = out.get_paths()[0];

        if (path.size() != 0) {
            path_len = path_length(path, tb.distance(test_label));
            std::cout << std::endl << "Path length: " << path_len << std::endl << std::endl;
        }
        else {
            std::cout << std::endl << "Path length: NaN" << std::endl << std::endl;
        }

        if (snapshot) {
            display_snapshots(tb.getws(test_label), path);
        }
        // plot graph
        else {
            if (test_label == "A" || test_label == "B" || test_label == "D" || test_label == "DT") {
                plot_graph(out, tb.test_colision(test_label), tb.add_obstacle_edges(test_label), tb.joint_limits(test_label), one_by_one, write_to_file);
            }
            else{
                plot_3d(out, tb.test_colision(test_label), tb.add_obstacle_edges(test_label), tb.joint_limits(test_label), one_by_one, write_to_file);
            }
        }
    }

    // if simulation mode is active run algorithm num_iter times and average the execution times
    else if (mode == "-sim") {
        // number of repeats
        size_t num_iter = 20;

        auto begin = std::chrono::high_resolution_clock::now();
        output out;

        for (size_t iter = 0; iter < num_iter; ++iter) {
            try {
                test_sq = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                    tb.get_sample(test_label), tb.distance(test_label), tb.num_samples(test_label), tb.stepsize(test_label), tb.radius(test_label));

                out = test_sq.run_test(algorithm);
            }
            catch (std::logic_error err) {
                std::cout << err.what() << std::endl;
                throw;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::cout << std::endl << "Average elapsed time over " << num_iter << " executions: " << elapsed.count() / num_iter * 1e-9 << " seconds" << std::endl << std::endl;
    }

    // simulate various parameter values
    else if (mode == "-param") {
        // number of iterations
        const size_t num_iter = 27;
        const size_t num_repeats = 5;
        std::vector <output> out(num_iter);
        std::vector <test> tests(num_iter);
        std::vector <double> elapsed_time(num_iter);
        std::vector <double> avg_length(num_iter);
        std::vector <bool> path_failed(num_iter);

        double radii[] = {tb.radius(test_label) / 2, tb.radius(test_label), tb.radius(test_label) * 2};
        double stepsizes[] = {tb.stepsize(test_label) / 2, tb.stepsize(test_label), tb.stepsize(test_label) * 2};
        unsigned int nums_samples[] = {tb.num_samples(test_label) / 2, tb.num_samples(test_label), tb.num_samples(test_label) * 2};

        try {
            for (size_t iter = 0; iter < num_iter; ++iter) {
                auto begin = std::chrono::high_resolution_clock::now();
                double radius, stepsize;
                unsigned int num_samples;
                
                std::cout << std::endl;
                std::cout << "*************************" << std::endl;
                std::cout << "Iteration " << iter+1 << std::endl;
                std::cout << "*************************" << std::endl;

                radius = radii[iter % 3];
                stepsize = stepsizes[(iter / 3) % 3];
                num_samples = nums_samples[iter / 9];

                tests[iter] = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                    tb.get_sample(test_label), tb.distance(test_label), num_samples, stepsize, radius);

                // calculate average path length
                for (size_t repeat = 0; repeat < num_repeats; ++repeat) {
                    // get path
                    out[iter] = tests[iter].run_test(algorithm);
                    auto path = out[iter].get_paths()[0];
                    double path_len = path_length(path, tb.distance(test_label));

                    // test existance of path
                    if (path.size() > 0) {
                        avg_length[iter] += path_len;
                    }
                    else {
                        path_failed[iter] = true;
                    }
                }

                avg_length[iter] /= num_repeats;

                auto end = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                elapsed_time[iter] = elapsed.count() / num_repeats * 1e-9;
            }

            // output list of parameter combinations and relevant diagnostics to file
            std::ofstream ofs("diagnostics" + algorithm + ".txt");
            ofs << "Diagnostic information" << std::endl;
            ofs << "Path length and runtime averaged over 5 executions" << std::endl;
            ofs << "Connected if path was found in all 5 executions" << std::endl << std::endl;

            for (size_t iter = 0; iter < num_iter; ++iter) {
                ofs << "test" << test_label << " (num_samples = " << tests[iter].get_num_samples();
                ofs << ", stepsize = " << tests[iter].get_stepsize();
                ofs << ", radius = " << tests[iter].get_radius() << ") : ";;

                if (!path_failed[iter]) {
                    ofs << "is connected (path_length = " << avg_length[iter] << ", elapsed_time = " << elapsed_time[iter] << ")" << std::endl;
                }
                else {
                    ofs << "is not connected" << std::endl;
                }
            }
        }
        catch (std::logic_error err) {
            std::cout << err.what() << std::endl;
            throw;
        }
    }

    return 0;
}