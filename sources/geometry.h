#pragma once

#include <vector>
#include <ostream>
#include <functional>
#include <iostream>
#include <cmath>

// class implementing point
template <typename real>
class point {
protected:
    std::vector <real> components;

public:
    point () {}

    point (const size_t dim) {
        (this->components).resize(dim);
    }

    point (const std::vector <real>& comp) : components(comp) {}

    point (const point <real>& P) : components(P.components) {}

    std::vector <real> get_components () const {
        return this->components;
    }

    size_t get_dimension () const {
        return (this->components).size();
    }

    real& operator[] (size_t index) {
        return this->components[index];
    }

    const real& operator[] (size_t index) const {
        return this->components[index];
    }

    point <real> operator+ (const point <real>& other) const {
        if (this->get_dimension() != other.get_dimension()) {
            throw std::length_error("Dimensions of points differ");
        }

        point <real> diff(this->get_dimension());

        for (size_t dim = 0; dim < this->get_dimension(); ++dim) {
            diff[dim] = this->operator[](dim) + other[dim];
        }
        
        return diff;
    }

    point <real> operator- (const point <real>& other) const {
        if (this->get_dimension() != other.get_dimension()) {
            throw std::length_error("Dimensions of points differ");
        }

        point <real> diff(this->get_dimension());

        for (size_t dim = 0; dim < this->get_dimension(); ++dim) {
            diff[dim] = this->operator[](dim) - other[dim];
        }
        
        return diff;
    }

    point <real> operator* (real scaler) const {
        point <real> scaled(this->get_dimension());

        for(size_t dim = 0; dim < this->get_dimension(); ++dim) {
            scaled[dim] = this->operator[](dim) * scaler;
        }

        return scaled;
    }   

    // implementing with safe floating point comparisson

    bool operator== (const point <real>& other) {
        const real epsilon = 1e-6;

        for (size_t dim = 0; dim < components.size(); ++dim) {
            real dif = components[dim] - other[dim];

            // difference more than epsilon
            if (dif < -epsilon || dif > epsilon) {
                return false;
            }
        }

        return true;
    }

    bool operator!= (const point <real>& other) {
        return !((*this) == other);
    }
};


template <typename real>
class point2d : public point <real> {
public:
    point2d () {}

    point2d (const point2d<real>& P) {
        point <real>::components = P.get_components();
    }

    point2d (const point<real>& P) {
        point <real>::components = P.get_components();
    }

    point2d (const double x, const double y) {
        point<real>::components.resize(2);
        point<real>::components[0] = x;
        point<real>::components[1] = y;
    }

    real getx () const {
        return point<real>::components[0];
    }

    real gety () const {
        return point<real>::components[1];
    }
};

template <typename real>
class point3d : public point <real> {
public:
    point3d () {}

    point3d (const point3d<real>& P) {
        point <real>::components = P.get_components();
    }

    point3d (const point<real>& P) {
        point <real>::components = P.get_components();
    }

    point3d (const double x, const double y, const double z) {
        point<real>::components.resize(3);
        point<real>::components[0] = x;
        point<real>::components[1] = y;
        point<real>::components[2] = z;
    }
    
    real getx () const {
        return point<real>::components[0];
    }

    real gety () const {
        return point<real>::components[1];
    }

    real getz () const {
        return point<real>::components[2];
    }
};

// axis-parallel box in N-d space
template <typename real>
class box {
private:
    std::vector <std::pair <real, real>> box_limits;

public:
    box () {}

    box (const box& B) : box_limits(B.box_limits) {}

    box (const std::vector <std::pair <real, real>>& _box_limits) : box_limits(_box_limits) {}

    // get meshgrid
    void meshgrid (
        std::vector <std::vector <real>>& X, 
        std::vector <std::vector <real>>& Y,
        std::vector <std::vector <real>>& Z,
        const real stepsize
    ) const {
        if (box_limits.size() != 3) {
            throw std::length_error("Dimension must be 3 for meshgrid");
        }

        X.clear();
        Y.clear();
        Z.clear();

        // iterate over all 
        for (real x = box_limits[0].first; x <= box_limits[0].second; x += stepsize) {
            std::vector <real> xr, yr, zr;

            for (real y = box_limits[1].first; y <= box_limits[1].second; y += stepsize) {
                for (real z = box_limits[2].first; z <= box_limits[2].second; z += stepsize) {
                    xr.push_back(x);
                    yr.push_back(y);
                    zr.push_back(z);
                }
            }
            
            X.push_back(xr);
            Y.push_back(yr);
            Z.push_back(zr);
        }
    }


};


// operator overload for outputing point<double>

std::ostream& operator << (std::ostream& os, const point <double> P);

// euclidean distance
// returns the distance squared, in order to generalize to abstact types
template <typename real>
real euclidean_distance (const point<real>& A, const point<real>& B) {
    if (A.get_dimension() != B.get_dimension()) {
        throw std::length_error("Dimensions of points differ");
    }

    size_t dimension = A.get_dimension();
    real total = 0;

    for (size_t index = 0; index < dimension; ++index) {
        total += (A[index] - B[index]) * (A[index] - B[index]);
    }

    return sqrt(total);
}

// weighed euclidean distance
// returns the distance squared, in order to generalize to abstact types
template <typename real>
real weighed_euclidean (
    const point <real>& A,
    const point <real>& B,
    const std::vector <real>& link_length
) {
    if (A.get_dimension() != B.get_dimension()) {
        throw std::length_error("Dimensions of points differ");
    }

    if (link_length.size() != A.get_dimension()) {
        throw std::length_error("Number of links differs from configuration space dimension");
    }

    size_t dimension = A.get_dimension();
    real total = 0;

    // calculate cumulative link lengths
    std::vector <real> cum_weight(dimension);
    cum_weight[dimension - 1] = link_length[dimension - 1];

    // mind that index has to be signed for the loop to work
    for (int index = dimension-2; index >= 0; --index) {
        cum_weight[index] = cum_weight[index + 1] + link_length[index];
    }

    // adjust metric
    for (size_t index = 0; index < dimension; ++index) {
        total += cum_weight[index] * (A[index] - B[index]) * (A[index] - B[index]);
    }

    return total;
}

// bind weighed euclidean to link lengths
template <typename real>
std::function <real(const point <real>&, const point <real>&)> bind_weighed_euclidean (
    const std::vector <real>& link_length
) {
    using namespace std::placeholders;
    return std::bind(weighed_euclidean <real>, _3, link_length);
}


point <double> walk (
    const point <double>& nearest, 
    const point <double>& sample, 
    const double scaler
);

bool path_clear(
    const point <double> A,
    const point <double> B,
    const std::function <bool(point<double>)>& collision_check,
    const double stepsize,
    const double epsilon = 1e-6
);

double path_length (
    std::vector <point <double>> path,
    const std::function <double(const point <double>&, const point <double>&)> distance
);


// test if point lies on line
template <typename real>
bool on_line (
    const point <real> endpoint,
    const std::pair <point <real>, point <real>>& line,
    const real epsilon = 1e-6
) {
    // if and only if dist(A,B) + dist(B,C) = dist(A,C)
    return abs(euclidean_distance(endpoint, line.first) + euclidean_distance(endpoint, line.second) - euclidean_distance(line.first, line.second)) < epsilon;
}

/*
Test if lines intersect

(y2-y1) / (x2-x1) = (y-y1) / (x-x1)
(y2-y1) (x-x1) = (y-y1) (x2-x1)
(y-y1) (x2-x1) + (y1-y2) (x-x1) = 0
(x2-x1) y + (y1-y2) x = (x2-x1) y1 + (y1-y2) x1
Ax + By = Ax1 + By1
Ax + By = C

Solve system
A1 x + B1 y = C1
A2 x + B2 y = C2
*/

template <typename real>
bool lines_intersect2d (
    const std::pair <point2d <real>, point2d <real>>& linea,
    const std::pair <point2d <real>, point2d <real>>& lineb,
    const real epsilon = 1e-6
) {
    real B1 = linea.second.getx() - linea.first.getx();
    real A1 = linea.first.gety() - linea.second.gety();
    real B2 = lineb.second.getx() - lineb.first.getx();
    real A2 = lineb.first.gety() - lineb.second.gety();
    real C1 = A1 * linea.first.getx() + B1 * linea.first.gety();
    real C2 = A2 * lineb.first.getx() + B2 * lineb.first.gety();

    // calculate determinant of system matrix
    real det = A1 * B2 - A2 * B1;
    real x, y;

    // std::cout << "linea = " << linea.first << "," << linea.second << std::endl;
    // std::cout << "lineb = " << lineb.first << "," << lineb.second << std::endl;

    // if determinant is nonzero we have unique intersection point
    if (abs(det) > epsilon) {
        // if point lies on both of the line segments
        real x = (B2 * C1 - B1 * C2) / det;
        real y = (A1 * C2 - A2 * C1) / det;
        // std::cout << "intersection = " << "(" << x << "," << y << ")" << std::endl; 

        bool betweenxA = (x <= linea.first.getx() && x >= linea.second.getx()) ||  (x <= linea.second.getx() && x >= linea.first.getx());
        bool betweenyA = (y <= linea.first.gety() && y >= linea.second.gety()) ||  (y <= linea.second.gety() && y >= linea.first.gety());
        bool betweenxB = (x <= lineb.first.getx() && x >= lineb.second.getx()) ||  (x <= lineb.second.getx() && x >= lineb.first.getx());
        bool betweenyB = (y <= lineb.first.gety() && y >= lineb.second.gety()) ||  (y <= lineb.second.gety() && y >= lineb.first.gety());
        return (betweenxA && betweenyA && betweenxB && betweenyB);
    }

    // otherwise lines are parallel or overlapping

    // test if overlapping
    if (on_line <real>(linea.first, lineb) || on_line <real>(linea.second, lineb)) {
        return true;
    }

    // otherwise parallel
    return false;
}


// arbitrary polygon in 2d space
template <typename real>
class polygon {
private:
    // given in counterclockwise order
    std::vector <std::pair <point2d <real>, point2d <real>>> edges;

public:
    polygon () {}

    polygon (const polygon& P) : edges(P.edges) {}

    polygon (const std::vector <std::pair <point2d <real>, point2d<real>>>& _edges) : edges(_edges) {}

    std::vector <std::pair <point2d <real>, point2d <real>>> get_edges () const {
        return edges;
    }
};


// class representing generic 3D arm manipulator
template <typename real>
class arm {
protected:
    point <real> base;
    std::vector <real> link_lengths;
    std::vector <std::pair <real,real>> joint_limits;

public:
    arm () {}

    arm (const arm <real>& robot) : base(robot.base), link_lengths(robot.link_lengths), joint_limits(robot.joint_limits) {}

    arm (
        const point <real> _base,
        const std::vector <real>& _link_lengths,
        const std::vector <std::pair <real,real>>& _joint_limits
    ) : base(_base), link_lengths(_link_lengths), joint_limits(_joint_limits) {}

    std::vector <std::pair <real,real>> get_joint_limits () const {
        return joint_limits;
    }

    std::vector <real> get_link_lengths () const {
        return link_lengths;
    }

    bool intersects (
        const box <real>& obstacle,
        const point <real>& configuration
    ) const {
        throw std::logic_error("nah");
    }

    virtual bool intersects (
        const polygon <real>& obstacle,
        const point <real>& configuration
    ) const {
        throw std::logic_error("Never should call this method in base class");
    }

    virtual std::vector <std::pair <point <real>, point <real>>> dir_kine (
        const point <real> configuration
    ) const {
        throw std::logic_error("nah");
    }
};


// helper function for converting pair of points to pair of 2D points
template <typename real>
std::pair <point2d <real>, point2d <real>> linkto2d (
    const std::pair <point <real>, point <real>> link
) {
    point2d <real> fi(link.first);
    point2d <real> se(link.second);
    std::pair <point2d <real>, point2d <real>> link2d{fi, se};
    return link2d;
}


// class representing 2D arm manipulator
template <typename real>
class arm2d : public arm <real> {
public:
    arm2d () {}

    arm2d (const arm2d <real>& robot) : arm<real>(robot.base, robot.link_lengths, robot.joint_limits) {}

    arm2d (
        const point2d <real> _base,
        const std::vector <real>& _link_lengths,
        const std::vector <std::pair <real,real>>& _joint_limits
    ) : arm<real>(_base, _link_lengths, _joint_limits) {}

    // check if 2D arm intersects polygon
    bool intersects (
        const polygon <real>& obstacle,
        const point <real>& configuration
    ) const override {
        auto links = dir_kine(configuration);

        // check intersection of links with obstacles
        for (const auto& link : links) {
            for (const auto& edge : obstacle.get_edges()) {
                std::pair <point2d <real>, point2d <real>> link2d = linkto2d(link);
                if (lines_intersect2d(link2d, edge)) {
                    return true;
                }
            }
        }

        // check intersection of links amongst themselves
        for (auto it = links.begin(); it != links.end(); ++it) {
            if (std::next(it, 1) != links.end()) {
                for (auto jt = std::next(it, 2); jt != links.end(); ++jt) {
                    std::pair <point2d <real>, point2d <real>> A = linkto2d(*it);
                    std::pair <point2d <real>, point2d <real>> B = linkto2d(*jt);
                    if (lines_intersect2d(A, B)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    // solve direct kinematics for planar arm
    std::vector <std::pair <point <real>, point <real>>> dir_kine (
        const point <real> configuration
    ) const override {

        if (configuration.get_dimension() != arm<real>::link_lengths.size()) {
            throw std::length_error("Dimension of configuration and number of joints in planar arm differ");
        }

        for (size_t dim = 0; dim < configuration.get_dimension(); ++dim) {
            if (configuration[dim] < arm<real>::joint_limits[dim].first || configuration[dim] > arm<real>::joint_limits[dim].second) {
                throw std::domain_error("Configuration coordinates lie outside joint limits");
            }
        }   

        // calculate links, based on configuration vector
        std::vector <std::pair <point <real>, point <real>>> links;
        point <real> current(arm<real>::base);
        real total_angle = 0;

        for (size_t dim = 0; dim < configuration.get_dimension(); ++dim) {
            total_angle += configuration[dim];
            point <real> next = current + point2d <real>(cos(total_angle), sin(total_angle)) * arm<real>::link_lengths[dim];
            // std::cout << "link(" << total_angle * 180 / M_PI << ") = " << "(" << current << "," << next << ")" << std::endl;
            links.push_back({current, next});
            current = next;
        }

        return links;
    } 
};

// class representing workspace
template <typename real>
class workspace {
protected:
    arm <real>* robot;

public:
    workspace () {}

    workspace (const workspace <real>& W) {
        robot = new arm <real>(*W.robot);
    }

    workspace (
        const arm <real>& _robot
    ) {
        robot = new arm <real>(_robot);
    }

    workspace (
        arm <real>* const _robot
    ) {
        robot = _robot;
    }

    // return link lengths
    std::vector <real> get_link_lengths () const {
        return robot->get_link_lengths();
    }

    arm <real>& get_robot () const {
        return *robot;
    }

    ~workspace () {
        delete robot;
    }
};


// define 2D workspace
template <typename real>
class workspace2d : public workspace <real> {
private:
    std::vector <polygon<real>> obstacles;

public:
    workspace2d () {}

    workspace2d (const workspace2d <real>& ws) : obstacles(ws.obstacles), workspace <real>(*ws.robot) {}

    workspace2d (
        const std::vector <polygon<real>> _obstacles,
        const arm <real>& _robot
    ) : obstacles(_obstacles), workspace <real>(_robot) {}

    workspace2d (
        const std::vector <polygon<real>> _obstacles,
        arm <real>* const _robot
    ) : obstacles(_obstacles), workspace <real>(_robot) {}

    // test if specific configuration colifs
    bool collides (const point <real> configuration) {
        for (const auto& obstacle : obstacles) {
            if (workspace<real>::robot->intersects(obstacle, configuration)) {
                return true;
            }
        }
        return false;
    }

    const std::vector <polygon<real>> get_obstacles () const {
        return obstacles;
    }
};


// define 3D workspace
template <typename real>
class workspace3d : public workspace <real> {
    
};