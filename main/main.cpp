#include "geometry.h"
#include "planning.h"
#include "rrt.h"
#include "prm.h"
#include "fmt.h"
#include "plot.h"
#include <functional>
#include <random>
#include <sstream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>

class test {
private:
    std::vector <point<double>> starts, goals;
    std::vector <std::pair <double,double>> joint_limits;
    unsigned int num_samples;
    double stepsize, radius;
    std::function<bool(point<double>)> test_collision;
    std::function<point<double>()> get_sample;
    std::function <double(const point <double>&, const point <double>&)> distance;
    double eta;

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
        const double _radius,
        const double _eta = 1
    ) : starts(_starts), goals(_goals), joint_limits(_joint_limits), num_samples(_num_samples), stepsize(_stepsize), radius(_radius),
    test_collision(_test_collision), get_sample(_get_sample), distance(_distance), eta(_eta) {

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
        const double _radius,
        const double _eta = 1
    ) : joint_limits(_joint_limits), num_samples(_num_samples), stepsize(_stepsize), radius(_radius),
    test_collision(_test_collision), get_sample(_get_sample), distance(_distance), eta(_eta) {

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
            out = rrt(start, goal, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else if (algorithm == "FMT") {
            out = fmt(start, goal, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else if (algorithm == "PRM") {
            out = prm(starts, goals, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, radius);
        }
        else if (algorithm == "FMT*") {
            out = fmtstar(start, goal, joint_limits, test_collision, get_sample, distance, num_samples, stepsize, eta);
        }
        else {
            throw std::domain_error("Algorithm label has to be either RRT, PRM, FMT or FMT*");
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
    constexpr static double radiusC = 0.1; // same as old

    constexpr static unsigned int num_samplesD = 300; // since there are more collisions
    constexpr static double stepsizeD = 3e-2; // around PI * old_value
    constexpr static double radiusD = 1.7; // since weighed euclidean, have to scale

    constexpr static unsigned int num_samplesE = 500;
    constexpr static double stepsizeE = 3e-2;
    constexpr static double radiusE = 4;

    constexpr static unsigned int num_samplesF = 20000;
    constexpr static double stepsizeF = 0.1;
    constexpr static double radiusF = 20;

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
        std::vector <point <double>>& previous
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
        return get_sample_util(joint_limitsA());
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
        std::vector <point <double>>& previous
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
        return get_sample_util(joint_limitsB());
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

    static box <double> box_downleft() {
        std::vector <std::pair <double, double>> lims{
            {0, 0.4}, {0, 0.4}, {0, 1}
        };
        return box <double>(lims);
    }

    static box <double> box_upright() {
        std::vector <std::pair <double, double>> lims{
            {0.6, 1}, {0.6, 1}, {0, 1}
        };
        return box <double>(lims);
    }

    static box <double> box_downright() {
        std::vector <std::pair <double, double>> lims{
            {0.6, 1}, {0, 0.4}, {0, 1}
        };
        return box <double>(lims);
    }

    static box <double> box_upleft() {
        std::vector <std::pair <double, double>> lims{
            {0, 0.4}, {0.6, 1}, {0, 1}
        };
        return box <double>(lims);
    }

    static workspace3d <double> workspaceC () {
        std::vector <box <double>> obstacles;
        obstacles.push_back(box_upleft());
        obstacles.push_back(box_upright());
        obstacles.push_back(box_downleft());
        obstacles.push_back(box_downright());

        arm <double>* planar = nullptr;
        workspace3d <double> ws(obstacles, planar);
        return ws;
    }

    static point <double> get_sampleC () {
        return get_sample_util(joint_limitsC());
    }

    static size_t add_obstacle_edgesC (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous 
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

    static workspace2d <double> workspaceD (const std::string label) {
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

        arm <double>* planar = new arm2d <double>(base, link_lengths, lims);
        workspace2d <double> ws(obstacles, planar);
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

    static std::function <double(const point <double>&, const point <double>&)> distanceD () {
        auto linklen = workspaceD("D").get_link_lengths();
        auto dlam = [linklen](const point <double>& A, const point <double>& B) {
            return weighed_euclidean <double>(A, B, linklen);
        };
        return dlam;
    }

    static std::function <double(const point <double>&, const point <double>&)> distanceDT () {
        auto linklen = workspaceD("DT").get_link_lengths();
        auto dlam = [linklen](const point <double>& A, const point <double>& B) {
            return weighed_euclidean <double>(A, B, linklen);
        };
        return dlam;
    }

    static std::function<bool(point<double>)> test_collisionD () {
        auto ws = workspaceD("D");
        auto clam = [ws](const point <double> config) {
            return ws.collides(config);
        };
        return clam;
    }

    static std::function<bool(point<double>)> test_collisionDT () {
        auto ws = workspaceD("DT");
        auto clam = [ws](const point <double> config) {
            return ws.collides(config);
        };
        return clam;
    }

    static point <double> get_sampleD () {
        return get_sample_util(joint_limitsD());
    }

    /***************************
    TEST E (antropomorphic arm)
    ***************************/

    static box <double> box_baseE () {
        // box as base (table, ground or whatever)
        std::vector <std::pair <double, double>> lims {
            {-15, 15}, {-15, 15}, {-15, 0}
        };
        return box <double>(lims);
    }

    static box <double> floating_boxE () {
        // floating box which will need to be avoided
        std::vector <std::pair <double, double>> lims {
            {-2, 2}, {-15, 15}, {6, 15}
        };
        return box <double>(lims);
    }

    static point3d <double> baseE () {
        point3d <double> base(0, 0, 1e-4);
        return base;
    }

    static std::vector <double> link_lengthsE () {
        std::vector <double> link_lengths{5, 5, 6};
        return link_lengths;
    }

    static std::vector <std::pair <double,double>> joint_limitsE () {
        std::vector <std::pair <double,double>> lims {{-M_PI, M_PI}, {-M_PI, M_PI}, {-M_PI, M_PI}};
        return lims;
    }

    static workspace3d <double> workspaceE () {
        std::vector <box <double>> obstacles;
        obstacles.push_back(box_baseE());
        obstacles.push_back(floating_boxE());

        point3d <double> base(baseE());
        std::vector <double> link_lengths{link_lengthsE()};
        std::vector <std::pair <double,double>> lims(joint_limitsE());

        arm <double>* antro = new antropomorphic_arm <double>(base, link_lengths, lims);
        workspace3d <double> ws(obstacles, antro);
        return ws;
    }

    static point <double> startE () {
        return point <double>(std::vector<double>{0, 0.2, 1.7});
    }

    static point <double> goalE () {
        return point <double>(std::vector<double>{M_PI, 0.2, 1.7});
    }

    static std::vector<point <double>> startsE () {
        std::vector <point <double>> starts{startE()};
        return starts;
    }

    static std::vector <point <double>> goalsE () {
        std::vector <point <double>> goals{goalE()};
        return goals;
    }

    static std::function <double(const point <double>&, const point <double>&)> distanceE () {
        auto linklen = workspaceE().get_link_lengths();
        auto dlam = [linklen](const point <double>& A, const point <double>& B) {
            return weighed_euclidean <double>(A, B, linklen);
        };
        return dlam;
    }

    static std::function<bool(point<double>)>test_collisionE () {
        auto ws = workspaceE();
        auto clam = [ws](const point <double> config) {
            return ws.collides(config);
        };
        return clam;
    }

    static point <double> get_sampleE () {
        return get_sample_util(joint_limitsE());
    }

    static size_t add_obstacle_edgesE (
        std::vector <point <double>>& current, 
        std::vector <point <double>>& previous 
    ) {
        // floating box

        previous.push_back(point2d <double>(-2, -2));
        current.push_back(point2d <double>(-2, 2));

        previous.push_back(point2d <double>(-2, 2));
        current.push_back(point2d <double>(2, 2));

        previous.push_back(point2d <double>(2, 2));
        current.push_back(point2d <double>(2, -2));
        
        previous.push_back(point2d <double>(2, -2));
        current.push_back(point2d <double>(-2, -2));

        // base

        previous.push_back(point2d <double>(-15, -15));
        current.push_back(point2d <double>(-15, 15));

        previous.push_back(point2d <double>(-15, 15));
        current.push_back(point2d <double>(15, 15));

        previous.push_back(point2d <double>(15, 15));
        current.push_back(point2d <double>(15, -15));
        
        previous.push_back(point2d <double>(15, -15));
        current.push_back(point2d <double>(-15, -15));

        return 8;
    }

    /***************************
    TEST F (planar arm)
    ***************************/

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesF_A () {
        // define points of obstacle
        const double hi = 8, lo = 6.5;
        point2d <double> A(-hi, 2);
        point2d <double> B(-hi, -2);
        point2d <double> C(-lo, -2);
        point2d <double> D(-lo, 2);
        point2d <double> E(-hi, 2);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesF_B () {
        // define points of obstacle
        const double hi = 2, lo = 0;
        point2d <double> A(-hi, 2);
        point2d <double> B(-hi, -2);
        point2d <double> C(-lo, -2);
        point2d <double> D(-lo, 2);
        point2d <double> E(-hi, 2);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesF_C () {
        // define points of obstacle
        const double hi = 2, lo = 0;
        point2d <double> A(hi, 2);
        point2d <double> B(hi, -2);
        point2d <double> C(lo, -2);
        point2d <double> D(lo, 2);
        point2d <double> E(hi, 2);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static std::vector <std::pair <point2d<double>, point2d<double>>> edgesF_D () {
        // define points of obstacle
        const double hi = 8, lo = 6.5;
        point2d <double> A(hi, 2);
        point2d <double> B(hi, -2);
        point2d <double> C(lo, -2);
        point2d <double> D(lo, 2);
        point2d <double> E(hi, 2);
        std::vector <point2d <double>> points{A, B, C, D, E};
        return edges_from_points(points);
    }

    static point2d <double> baseF () {
        point2d <double> base(0, 4.5);
        return base;
    }

    static std::vector <double> link_lengthsF () {
        std::vector <double> link_lengths{4, 4, 4, 4, 4, 4};
        return link_lengths;
    }

    static std::vector <std::pair <double,double>> joint_limitsF () {
        std::vector <std::pair <double,double>> lims {
            {-2*M_PI, 2*M_PI}, {-2*M_PI, 2*M_PI}, {-2*M_PI, 2*M_PI},
            {-2*M_PI, 2*M_PI}, {-2*M_PI, 2*M_PI}, {-2*M_PI, 2*M_PI}
        };
        return lims;
    }

    static workspace2d <double> workspaceF () {
        // form polygons from points
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesA(edgesF_A());
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesB(edgesF_B());
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesC(edgesF_C());
        std::vector <std::pair <point2d<double>, point2d<double>>> edgesD(edgesF_D());
        std::vector <polygon <double>> obstacles;

        obstacles.push_back(polygon <double>(edgesA));
        obstacles.push_back(polygon <double>(edgesB));
        obstacles.push_back(polygon <double>(edgesC)); 
        obstacles.push_back(polygon <double>(edgesD));

        // create 2D planar arm object
        point2d <double> base(baseF());
        std::vector <double> link_lengths{link_lengthsF()};
        std::vector <std::pair <double,double>> lims(joint_limitsF());

        arm <double>* planar = new arm2d <double>(base, link_lengths, lims);
        workspace2d <double> ws(obstacles, planar);
        return ws;
    }

    static point <double> startF () {
        return point <double>(std::vector<double>{-M_PI_4, -M_PI_4, 0, 0, 0, -M_PI_4});
    }

    static point <double> goalF () {
        return point <double>(std::vector<double>{5*M_PI_4, M_PI_4, 0, 0, 0, M_PI_4});
    }

    static std::vector<point <double>> startsF () {
        std::vector <point <double>> starts{startF()};
        return starts;
    }

    static std::vector <point <double>> goalsF () {
        std::vector <point <double>> goals{goalF()};
        return goals;
    }

    static std::function <double(const point <double>&, const point <double>&)> distanceF () {
        // create lambda and bind link lengths
        auto linklen = workspaceF().get_link_lengths();
        auto dlam = [linklen](const point <double>& A, const point <double>& B) {
            return weighed_euclidean <double>(A, B, linklen);
        };
        return dlam;
    }

    static std::function<bool(point<double>)> test_collisionF () {
        auto ws = workspaceF();
        auto clam = [ws](const point <double> config) {
            return ws.collides(config);
        };
        return clam;
    }

    static point <double> get_sampleF () {
        return get_sample_util(joint_limitsF());
    }

    /**********************
    HELPER FUNCTIONS
    ***********************/

    // functions returning test info with test label as argument
    static point <double> start (const std::string label) {
        if (label == "A") return startA();
        if (label == "B") return startB();
        if (label == "C") return startC();
        if (label == "D" || label == "DT") return startD();
        if (label == "E") return startE();
        if (label == "F") return startF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static point <double> goal (const std::string label) {
        if (label == "A") return goalA();
        if (label == "B") return goalB();
        if (label == "C") return goalC();
        if (label == "D") return goalD();
        if (label == "DT") return goalDT();
        if (label == "E") return goalE();
        if (label == "F") return goalF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::vector<point <double>> starts (const std::string label) {
        if (label == "A") return startsA();
        if (label == "B") return startsB();
        if (label == "C") return startsC();
        if (label == "D" || label == "DT") return startsD();
        if (label == "E") return startsE();
        if (label == "F") return startsF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::vector <point <double>> goals (const std::string label) {
        if (label == "A") return goalsA();
        if (label == "B") return goalsB();
        if (label == "C") return goalsC();
        if (label == "D") return goalsD();
        if (label == "DT") return goalsDT();
        if (label == "E") return goalsE();
        if (label == "F") return goalsF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::vector <std::pair<double,double>> joint_limits (const std::string label) {
        if (label == "A") return joint_limitsA();
        if (label == "B") return joint_limitsB();
        if (label == "C") return joint_limitsC();
        if (label == "D" || label == "DT") return joint_limitsD();
        if (label == "E") return joint_limitsE();
        if (label == "F") return joint_limitsF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static std::function <double(const point <double>&, const point <double>&)> distance (const std::string label) {
        if (label == "A" || label == "B" || label == "C") return euclidean_distance <double>;
        if (label == "D") return distanceD();
        if (label == "DT") return distanceDT();
        if (label == "E") return distanceE();
        if (label == "F") return distanceF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static unsigned int num_samples (const std::string label) {
        if (label == "A") return num_samplesA;
        if (label == "B") return num_samplesB;
        if (label == "C") return num_samplesC;
        if (label == "D" || label == "DT") return num_samplesD;
        if (label == "E") return num_samplesE;
        if (label == "F") return num_samplesF;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static double stepsize (const std::string label) {
        if (label == "A") return stepsizeA;
        if (label == "B") return stepsizeB;
        if (label == "C") return stepsizeC;
        if (label == "D" || label == "DT") return stepsizeD;
        if (label == "E") return stepsizeE;
        if (label == "F") return stepsizeF;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static double radius (const std::string label) {
        if (label == "A") return radiusA;
        if (label == "B") return radiusB;
        if (label == "C") return radiusC;
        if (label == "D" || label == "DT") return radiusD;
        if (label == "E") return radiusE;
        if (label == "F") return radiusF;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    std::function<bool(point<double>)> test_colision (const std::string label) {
        if (label == "A") return test_collisionA;
        if (label == "B") return test_collisionB;
        if (label == "C") return test_collisionC;
        if (label == "D") return test_collisionD();
        if (label == "DT") return test_collisionDT();
        if (label == "E") return test_collisionE();
        if (label == "F") return test_collisionF();

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    std::function<point<double>()> get_sample (const std::string label) {
        if (label == "A") return get_sampleA;
        if (label == "B") return get_sampleB;
        if (label == "C") return get_sampleC;
        if (label == "D" || label == "DT") return get_sampleD;
        if (label == "E") return get_sampleE;
        if (label == "F") return get_sampleF;

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    plot_object add_obstacle_edges (
        const std::string label,
        const bool obs = true
    ) {
        if (label == "A") return plot_object(add_obstacle_edgesA);
        if (label == "B") return plot_object(add_obstacle_edgesB);
        if (obs && (label == "D" || label == "DT")) return plot_object(joint_limitsD(), test_collisionD());
        if (!obs && (label == "D" || label == "DT")) return plot_object(joint_limitsD(), test_collisionD(), 1);  

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    plot_object obstacle_edges_3d (
        const std::string label,
        const bool obs = true
    ) {
        if (label == "C") return plot_object(workspaceC());
        if (obs && label == "E") return plot_object(joint_limitsE(), test_collisionE(), 10);
        if (!obs && label == "E") return plot_object(joint_limitsE(), test_collisionE(), 1);

        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static workspace2d <double> getws (std::string label) {
        if (label == "D" || label == "DT") return workspaceD(label);
        if (label == "F") return workspaceF();
        
        // invalid test label
        throw std::domain_error("Invalid test label");
    }

    static workspace3d <double> getws3 (std::string label) {
        if (label == "C") return workspaceC();
        if (label == "E") return workspaceE();

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
    static point <double> get_sample_util (const std::vector <std::pair <double,double>>& joint_limits) {
        size_t dimension = joint_limits.size();
        std::vector <double> components(dimension);
        std::random_device r;

        for (size_t dim = 0; dim < dimension; ++dim) {
            double lower = joint_limits[dim].first;
            double upper = joint_limits[dim].second;
            std::uniform_real_distribution<double> unif(lower, upper);
            std::default_random_engine re(r());
            components[dim] = unif(re);
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

int main (int argc, char *argv[]) {
    test_battery tb;
    test test_sq, test_mq;

    std::string algorithm = "RRT";
    std::string test_label = "F";
    bool one_by_one = false;
    bool write_to_file = false;
    bool snapshot = false;
    std::string mode = "-normal";
    bool show_path = false;
    bool obs = true;
    double eta = 1;

    for (size_t i = 0; i < argc; ++i) {
        std::string flag(argv[i]);
        if (flag == "RRT" || flag == "PRM" || flag == "FMT" || flag == "FMT*") {
            algorithm = flag;
        }
        if (flag == "A" || flag == "B" || flag == "C" || flag == "D" || flag == "DT" || flag == "E" || flag == "F") {
            test_label = flag;
        }
        if (flag == "point2D-A") test_label = "A";
        if (flag == "point2D-B") test_label = "B";
        if (flag == "point3D") test_label = "C";
        if (flag == "arm2D") test_label = "D";
        if (flag == "trivial2D") test_label = "DT";
        if (flag == "antro") test_label = "E";
        if (flag == "hard") test_label = "F";
        if (flag == "-seq") {
            one_by_one = true;
        }
        if (flag == "-file") {
            write_to_file = true;
        }
        if (flag == "-timesim" || flag == "-paramsim" || flag == "-energysim" || flag == "-reachsim" || flag == "-normal") {
            mode = flag;
        }
        if (flag == "-snapshot") {
            snapshot = true;
        }
        if (flag == "-path") {
            show_path = true;
        }
        if (flag == "-nobs") {
            obs = false;
        }
        if (flag.substr(0, 4) == "eta=") {
            std::string num = flag.substr(4);
            eta = stod(num);
        }
    }

    // run algorithm normally and display results

    if (mode == "-normal") {
        auto begin = std::chrono::high_resolution_clock::now();
        output out;

        try {
            test_sq = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                tb.get_sample(test_label), tb.distance(test_label), tb.num_samples(test_label), tb.stepsize(test_label), tb.radius(test_label), eta);

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
            if (test_label == "D" || test_label == "DT" || test_label == "F") {
                display_snapshots(tb.getws(test_label), path);
            }
            else if (test_label == "E") {
                display_snapshots3d(tb.getws3(test_label), path);
            }
            else {
                std::cout << "Cannot take snapshots for test " << test_label << std::endl;
            }
        }
        // plot graph
        else {
            if (test_label == "A" || test_label == "B" || test_label == "D" || test_label == "DT") {
                plot_graph(out, tb.test_colision(test_label), tb.add_obstacle_edges(test_label, obs), tb.joint_limits(test_label), one_by_one, write_to_file);
            }
            else if (test_label == "C" || test_label == "E") {
                plot3d(out, tb.obstacle_edges_3d(test_label, obs), tb.joint_limits(test_label), show_path, write_to_file);
            }
            else {
                std::cout << "Cannot plot configuration space for test " << test_label << std::endl;
            }
        }
    }

    // if simulation mode is active run algorithm num_iter times and average the execution times
    else if (mode == "-timesim") {
        // number of repeats
        size_t num_iter = 20;

        auto begin = std::chrono::high_resolution_clock::now();
        output out;

        for (size_t iter = 0; iter < num_iter; ++iter) {
            try {
                test_sq = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                    tb.get_sample(test_label), tb.distance(test_label), tb.num_samples(test_label), tb.stepsize(test_label), tb.radius(test_label), eta);

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
    else if (mode == "-paramsim") {
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

                num_samples = nums_samples[iter % 3];
                stepsize = stepsizes[(iter / 3) % 3];
                radius = radii[iter / 9];

                tests[iter] = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                    tb.get_sample(test_label), tb.distance(test_label), num_samples, stepsize, radius, eta);

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

    else if (mode == "-energysim") {
        // number of repeats
        std::vector <double> samplecnt_scaler; 
        const double lower = 0.2;
        const double upper = 5;
        const double jump = 0.4;
        for (double step = lower; step <= upper; step += jump) {
            samplecnt_scaler.push_back(step);
        }
        std::vector <double> xs, ys;

        output out;
        const size_t num_iter = 3;

        for (size_t level = 0; level < samplecnt_scaler.size(); ++level) {
            size_t num_samples = tb.num_samples(test_label) * samplecnt_scaler[level];
            double energy = 0;
            unsigned int pathcnt = 0;

            for (size_t iter = 0; iter < num_iter; ++iter) {
                try {
                    test_sq = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                        tb.get_sample(test_label), tb.distance(test_label), num_samples, tb.stepsize(test_label), tb.radius(test_label), eta);

                    out = test_sq.run_test(algorithm);
                    auto path = out.get_paths()[0];
                    if (path.empty()) {
                        energy = NAN;
                    }
                    else {
                        energy += path_length(path, tb.distance(test_label));
                    }
                }
                catch (std::logic_error err) {
                    std::cout << err.what() << std::endl;
                    throw;
                }
            }
        
            xs.push_back(samplecnt_scaler[level]);
            ys.push_back(energy / num_iter);
        }

        plot_function(xs, ys, "number of samples", "energy", "Energy vs. Number of samples");
    }

    else if (mode == "-reachsim") {
        // number of repeats
        std::vector <double> samplecnt_scaler{0.01, 0.05, 0.1, 0.2, 0.5, 0.7, 1, 1.5, 2}; 
        std::vector <double> xs, ys;

        output out;
        const size_t num_iter = 10;

        for (size_t level = 0; level < samplecnt_scaler.size(); ++level) {
            size_t num_samples = tb.num_samples(test_label) * samplecnt_scaler[level];
            double pathcnt = 0;

            for (size_t iter = 0; iter < num_iter; ++iter) {
                try {
                    test_sq = test(tb.start(test_label), tb.goal(test_label), tb.joint_limits(test_label), tb.test_colision(test_label), 
                        tb.get_sample(test_label), tb.distance(test_label), num_samples, tb.stepsize(test_label), tb.radius(test_label), eta);

                    out = test_sq.run_test(algorithm);
                    auto path = out.get_paths()[0];
                    if (!path.empty()) {
                        pathcnt += 1;
                    }
                }
                catch (std::logic_error err) {
                    std::cout << err.what() << std::endl;
                    throw;
                }
            }
        
            xs.push_back(samplecnt_scaler[level]);
            ys.push_back(pathcnt / num_iter);
        }

        plot_function(xs, ys, "number of samples", "probability", "Probability of finding path vs. Number of samples");
    }

    return 0;
}