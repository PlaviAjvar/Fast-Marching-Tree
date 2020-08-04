#pragma once

#include <vector>
#include <ostream>
#include <functional>
#include <iostream>
#include <cmath>
#include <memory>

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

    real operator[] (size_t index) const {
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

    point <real> operator/ (real scaler) const {
        if (scaler == 0) {
            throw std::domain_error("Cannot divide by zero");
        }
        return this->operator* (1 / scaler);
    }

    // implementing with safe floating point comparisson

    bool operator== (const point <real>& other) const {
        if (components.size() != other.get_dimension()) {
            return false;
        }
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

    bool operator!= (const point <real>& other) const {
        return !((*this) == other);
    }

    // scalar product
    real operator* (const point <real>& other) const {
        if (components.size() != other.get_dimension()) {
            throw std::length_error("Vectors have to have same dimensions in dot product");
        }

        real total = 0;
        for (size_t dim = 0; dim < components.size(); ++dim) {
            total += components[dim] * other[dim];
        }
        return total;
    }

    real magnitude () const {
        real total = 0;
        for (size_t dim = 0; dim < components.size(); ++dim) {
            total += components[dim] * components[dim];
        }
        return sqrt(total);
    }

    point <real> normalize (const real epsilon = 1e-6) const {
        if (magnitude() < epsilon) return *this;
        return (*this) / magnitude();
    }

    // erase coordinate from point
    void erase (size_t dim) {
        components.erase(components.begin() + dim);
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

    point2d (const real x, const real y) {
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

    point3d (const std::vector <real>& _components) : point<real>(_components) {}
    
    real getx () const {
        return point<real>::components[0];
    }

    real gety () const {
        return point<real>::components[1];
    }

    real getz () const {
        return point<real>::components[2];
    }

    // cross product
    point3d <real> operator^ (const point3d <real>& other) const {
        std::vector <real> a(point <real>::components), b(other.components);
        std::vector <real> comp{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
        return point3d <real>(comp);
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

    // get edges which define box
    std::vector <std::pair <point <double>, point <double>>> get_edges () const {
        size_t dimension = box_limits.size();
        std::vector <point <real>> corners(1 << dimension);

        // first obtain all corners
        // bitmask all possible combinations
        for (size_t mask = 0; mask < (1 << dimension); ++mask) {
            std::vector <real> components(dimension);

            for (size_t dim = 0; dim < dimension; ++dim) {
                if (mask & (1 << dim)) {
                    components[dim] = box_limits[dim].second;
                }
                else {
                    components[dim] = box_limits[dim].first;
                }
            }

            corners[mask] = point <real>(components);
        }

        // obtain edges from corners
        std::vector <std::pair <point <double>, point <double>>> edges((1 << dimension) * dimension);

        for (size_t mask = 0; mask < (1 << dimension); ++mask) {
            for (size_t dim = 0; dim < dimension; ++dim) {
                // flip one bit in mask
                size_t nask = mask ^ (1 << dim);
                size_t idx = mask * dimension + dim;
                edges[idx] = {corners[mask], corners[nask]};
            }
        }

        return edges;
    }

    // get faces of box
    std::vector <std::pair <real, real>> get_box_limits () const {
        return box_limits;
    }
};

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


// operator overload for outputing point<double>

std::ostream& operator << (std::ostream& os, const point <double> P);

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

    return sqrt(total);
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
    const real epsilon = 1e-10
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

    // if determinant is nonzero we have unique intersection point
    if (abs(det) > epsilon) {
        // if point lies on both of the line segments
        real x = (B2 * C1 - B1 * C2) / det;
        real y = (A1 * C2 - A2 * C1) / det;
        point2d <real> inter(x, y);

        return (on_line<real>(inter, linea) && on_line<real>(inter, lineb));
    }

    // otherwise lines are overlapping (parallel is impossible because of zero distance between lines)
    if (on_line<real>(linea.first, lineb) || on_line<real>(linea.second, lineb) || on_line<real>(lineb.first, linea) || on_line<real>(lineb.second, linea)) {
        return true;
    }
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

    virtual bool intersects (
        const polygon <real>& obstacle,
        const point <real>& configuration
    ) const {
        throw std::logic_error("Never should call this method in base class");
    }

    virtual bool intersects (
        const box <real>& obstacle,
        const point <real>& configuration
    ) const {
        throw std::logic_error("Never should call this method in base class");
    }

    virtual bool intersects_self (
        const point <real>& configuration
    ) const {
        throw std::logic_error("Never should call this method in base class");
    }

    virtual std::vector <std::pair <point <real>, point <real>>> dir_kine (
        const point <real> configuration
    ) const {
        throw std::logic_error("Never should call this method in base class");
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

// helper function for converting pair of points to pair of 3D points
template <typename real>
std::pair <point2d <real>, point2d <real>> linkto3d (
    const std::pair <point <real>, point <real>> link
) {
    point3d <real> fi(link.first);
    point3d <real> se(link.second);
    std::pair <point3d <real>, point3d <real>> link3d{fi, se};
    return link3d;
}


// check if point is in box
template <typename real>
bool in_box_bounds (
    point <real> sample,
    const std::vector <std::pair <real, real>>& box_limits,
    const real epsilon = 1e-6
) {
    for (size_t dim = 0; dim < sample.get_dimension(); ++dim) {
        if (sample[dim] < box_limits[dim].first - epsilon || sample[dim] > box_limits[dim].second + epsilon) {
            return false;
        }
    }
    
    return true;
}


/* 
Check if 3d line intersercts 3d box

x = x0 + rx * t
y = y0 + ry * t
z = z0 + rz * t

p = p0 + r * t

-----------------------------------

Intersect line:
u = u0 + ru * t
with plane:
u = up

t = (up - u0) / ru
*/

template <typename real>
bool line_intersects_box (
    const std::pair <point <real>, point <real>>& link,
    const std::vector <std::pair <real, real>>& box_limits,
    const real epsilon = 1e-6
) {
    if (link.first.get_dimension() != box_limits.size() || link.second.get_dimension() != box_limits.size()) {
        throw std::length_error("Dimension of link differs from dimension of box");
    }

    // line is inside box
    if (in_box_bounds(link.first, box_limits) || in_box_bounds(link.second, box_limits)) {
        return true;
    }

    for (size_t dim = 0; dim < box_limits.size(); ++dim) {
        // get base poit and direction vector
        point <real> p0 = link.first;
        point <real> r = link.second - link.first;

        // check intersections with lower and upper side
        point <real> low_int, high_int;

        // standard case: direction vector is nonzero in this direction 

        if (abs(r[dim]) > epsilon) {
            real t_low = (box_limits[dim].first - p0[dim]) / r[dim];
            real t_high = (box_limits[dim].second - p0[dim]) / r[dim];

            // calculate intersection using thse parameters
            low_int = p0 + r * t_low;
            high_int = p0 + r * t_high;

            // check if these points are inside box limits, and that the intersection is not edge point
            if (t_low >= -epsilon && t_low <= 1+epsilon && in_box_bounds(low_int, box_limits)) {
                return true;
            }

            if (t_high >= -epsilon && t_high <= 1+epsilon && in_box_bounds(high_int, box_limits)) {
                return true;
            }
        }

        // otherwise, line is parallel with face of box
        // if line intersects face, then line will also intersect a perpendicular face
        // therefore there is no reason to check this case
    }
    
    return false;
}

template <typename real>
real magnitude (
    const point <real> point
) {
    return point.magnitude();
}


/* 
Find distance between lines in 3D space.
The shortast distance is along vector v perpendicular to both lines.

v ^ r = v ^ s = 0 ==> v ~ r ^ s

Then project vector u (between points on two lines), onto vector v to find shortest length.

distance = (u * v) / |v|

*/
template <typename real>
real line_distance (
    const point3d <real> r,
    const point3d <real> s,
    const point3d <real> u
) {
    point3d <real> v = r ^ s;
    
    // if r and s are parallel
    if (v == point3d <real>(0, 0, 0)) {
        // calculate distance as area of rhombus |u ^ r| over base of rhombus |r|
        return magnitude(u ^ r) / magnitude(r);
    }

    // otherwise they aren´t parallel
    // find projection of vector u
    // onto vector v, which is perpendicular to both r and s
    return (u * v) / v.magnitude();
}


// if colinear vectors point in same direction
template <typename real>
bool same_direction (
    const point <real> v,
    const point <real> u,
    const real epsilon = 1e-6
) {
    if (magnitude(v) < epsilon || magnitude(u) < epsilon) return false;
    return v.normalize() == u.normalize();
}


// get a scaler along direction vector
template <typename real>
real signed_scaler (
    const point3d <real>& vector,
    const point3d <real>& dir
) {
    real t = magnitude(vector) / magnitude(dir);

    if (!same_direction(vector, dir)) {
        t = -t;
    }

    return t;
}

/* 
Check if lines intersect in 3d

p = p0 + r * t
q = q0 + s * t

x = x0 + r * t1
x = x1 + s * t2

x1 - x0 = rx * t1 - sx * t2
y1 - y0 = ry * t1 - sy * t2
z1 - z0 = rz * t1 - sz * t2
*/

template <typename real>
bool lines_intersect3d(
    const std::pair <point3d <real>, point3d <real>> linea,
    const std::pair <point3d <real>, point3d <real>> lineb,
    const real epsilon = 1e-6
) {
    // get base poit and direction vector for both lines
    point3d <real> p0 = linea.first;
    point3d <real> r = linea.second - linea.first;
    point3d <real> q0 = lineb.first;
    point3d <real> s = lineb.second - lineb.first;

    // vector between lines
    point3d <real> u = q0 - p0;

    // check that distance between lines is zero
    if (line_distance(r, s, u) > epsilon) {
        return false;
    }

    // test if lines are overlapping (cannot be parallel because of zero distance)
    if (magnitude(r ^ s) < epsilon) {
        real t_fi = signed_scaler<real>(lineb.first - p0, r);
        real t_se = signed_scaler<real>(lineb.second - p0, r);

        // check if one of these fallw between 0 and 1, or if segment contains other segment
        return ((t_fi >= 0 && t_fi <= 1) || (t_se >= 0 && t_se <= 1) || (t_fi < 0 && t_se > 1) || (t_se < 0 && t_fi > 1));
    }

    // find distance scaler
    real scaler = signed_scaler<real>(s ^ u, s ^ r);

    // find new point
    point3d <real> inter = p0 + r * scaler;

    // check if new point lies between limits
    real t1 = signed_scaler<real>(inter - p0, r);
    real t2 = signed_scaler<real>(inter - q0, s);

    return (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1);
}


// class implementing matrices
template <typename real>
class matrix {
private:
    std::vector <std::vector <real>> mat;

public:
    matrix () {}

    matrix (const matrix <real>& other) : mat(other.mat) {}

    matrix (
        size_t rows,
        size_t cols
    ) {
        mat.resize(rows);
        for (auto& row : mat) {
            row.resize(cols);
        }
    }

    matrix (const std::vector <std::vector <real>>& _mat) : mat(_mat) {}

    // construct from point, give homogeneus coordinate form
    matrix (
        const point <real> vector
    ) {
        mat.resize(vector.get_dimension() + 1);

        for (size_t dim = 0; dim < vector.get_dimension(); ++dim) {
            mat[dim].resize(1);
            mat[dim][0] = vector[dim];
        }

        (mat.back()).resize(1);
        (mat.back())[0] = 1;
    }

    size_t rowcount () const {
        return mat.size();
    }

    size_t colcount () const {
        return mat[0].size();
    }

    std::vector <real>& operator[] (size_t ridx) {
        return mat[ridx];
    }

    std::vector <real> operator[] (size_t cidx) const {
        return mat[cidx];
    }

    matrix <real> operator+ (const matrix <real> other) const {
        if (rowcount() != other.rowcount() || colcount() != other.colcount()) {
            throw std::length_error("Matrix dimensions must be the same");
        }

        matrix <real> result(rowcount(), colcount());

        for (size_t ridx = 0; ridx < rowcount(); ++ridx) {
            for (size_t cidx = 0; cidx < colcount(); ++cidx) {
                result.mat[ridx][cidx] = mat[ridx][cidx] + other[ridx][cidx];
            }
        }

        return result;
    }

    matrix <real> operator- (const matrix <real> other) const {
        if (rowcount() != other.rowcount() || colcount() != other.colcount()) {
            throw std::length_error("Matrix dimensions must be the same");
        }

        matrix <real> result(rowcount(), colcount());

        for (size_t ridx = 0; ridx < rowcount(); ++ridx) {
            for (size_t cidx = 0; cidx < colcount(); ++cidx) {
                result.mat[ridx][cidx] = mat[ridx][cidx] - other[ridx][cidx];
            }
        }

        return result;
    }

    matrix <real> operator* (const matrix <real> other) const {
        if (colcount() != other.rowcount()) {
            throw std::length_error("Number of columns of first matrix must equal the number of rows of the second");
        }

        matrix <real> result(rowcount(), other.colcount());

        for (size_t ridx = 0; ridx < rowcount(); ++ridx) {
            for (size_t cidx = 0; cidx < other.colcount(); ++cidx) {
                for (size_t k = 0; k < colcount(); ++k) {
                    result[ridx][cidx] += mat[ridx][k] * other[k][cidx];
                }
            }
        }

        return result;
    }

    matrix <real> operator* (const real scaler) const {
        matrix <real> result(rowcount(), colcount());

        for (size_t ridx = 0; ridx < rowcount(); ++ridx) {
            for (size_t cidx = 0; cidx < colcount(); ++cidx) {
                result[ridx][cidx] = mat[ridx][cidx] * scaler;
            }
        }

        return result;
    }

    // get point from homogeneous coordinates
    point <real> to_point () const {
        if (mat.size() != 4 || mat[0].size() > 1 ) {
            throw std::length_error("Invalid form for extracting point");
        }

        std::vector <real> vec{mat[0][0], mat[1][0], mat[2][0]};
        return point <real>(vec);
    }

    friend std::ostream& operator<< (std::ostream& os, const matrix <real> matr) {
        for (const auto& row : matr.mat) {
            for (const auto& element : row) {
                os << element << " ";
            }
            os << std::endl;
        }
        return os;
    }
};

// class representing DH parameters
template <typename real>
class denavit_hartenberg {
private:
    // a, d, phi, alpha
    std::vector <real> param;
    std::vector <bool> var;

public:
    denavit_hartenberg () {
        param.resize(4);
        var.resize(4);
    }

    denavit_hartenberg (const denavit_hartenberg <real>& other) 
    : param(other.param), var(other.var) {}

    denavit_hartenberg (
        const real _a,
        const real _d,
        const real _phi,
        const real _alpha,
        const std::vector <bool>& _var
    ) : param{_a, _d, _phi, _alpha}, var(_var) {
        if (var.size() != 4) {
            throw std::length_error("Var has to have 4 parameters");
        }
    }

    denavit_hartenberg (
        const std::vector <real>& _param,
        std::string type
    ) {
        param = _param;
        var.resize(4);
        var[1] = (type == "prismatic");
        var[2] = (type == "revolute");
    }

    real getpar (size_t idx) {
        if (var[idx]) {
            throw std::logic_error("Parameter a is a variable");
        }
        return param[idx];
    }

    real geta () const {
        return getpar(0);
    }

    real getd() const {
        return getpar(1);
    }

    real getphi () const {
        return getpar(2);
    }

    real getalpha() const {
        return getpar(3);
    }

    // get homogeneus transformation matrix from dh parameters
    matrix <real> get_matrix (const std::vector <real>& par) const {
        size_t par_count = 0;
        std::vector <real> parinst(param);

        for (size_t i = 0; i < var.size(); ++i) {
            if (var[i]) {
                if (par_count >= par.size()) {
                    throw std::length_error("Too little paramaters supplied");
                }
                parinst[i] = par[par_count];
                par_count++;
            }
        }

        if (par_count < par.size()) {
            throw std::length_error("Too many parameters supplied");
        }

        real a = parinst[0];
        real d = parinst[1];
        real phi = parinst[2];
        real alpha = parinst[3];

        std::vector <std::vector <real>> mat{
            {cos(phi), -sin(phi)*cos(alpha), sin(phi)*sin(alpha), a*cos(phi)},
            {sin(phi), cos(phi)*cos(alpha), -cos(phi)*sin(alpha), a*sin(phi)},
            {0, sin(alpha), cos(alpha), d},
            {0, 0, 0, 1}
        };

        return matrix<real>(mat);
    }
};


// class representing antropomorphic arm
template <typename real>
class arm3d : public arm <real> {
protected:
    std::vector <denavit_hartenberg <real>> dh;

public:
    arm3d () {}

    arm3d (const arm3d <real>& robot) : arm <real>(robot.base, robot.link_lengths, robot.joint_limits), dh(robot.dh) {}

    // pass DH parameters to constructor
    arm3d (
        const point <real> _base,
        const std::vector <real>& _link_lengths,
        const std::vector <std::pair <real,real>>& _joint_limits,
        std::vector <denavit_hartenberg <real>> _dh
    ) : arm<real>(_base, _link_lengths, _joint_limits), dh(_dh) {}

    arm3d (
        const point <real> _base,
        const std::vector <real>& _link_lengths,
        const std::vector <std::pair <real,real>>& _joint_limits
    ) : arm<real>(_base, _link_lengths, _joint_limits) {}

    // check if 3D arm intersects a 3D box
    bool intersects (
        const box <real>& obstacle,
        const point <real>& configuration
    ) const override {
        std::vector <std::pair <point <real>, point <real>>> links = dir_kine(configuration);
        std::vector <std::pair <real, real>> box_limits = obstacle.get_box_limits();

        // chek intersection of links with obstacles
        for (const auto& link : links) {
            if (line_intersects_box(link, box_limits)) {
                return true;
            }
        }

        return false;
    }

    // check if 3D arm intersects itself
    bool intersects_self (
        const point <real>& configuration
    ) const override {
        std::vector <std::pair <point <real>, point <real>>> links = dir_kine(configuration);

        // check intersection of links amongst themselves
        for (auto it = links.begin(); it != links.end(); ++it) {
            if (std::next(it, 1) != links.end()) {
                for (auto jt = std::next(it, 2); jt != links.end(); ++jt) {
                    std::pair <point3d <real>, point3d <real>> A = linkto3d(*it);
                    std::pair <point3d <real>, point3d <real>> B = linkto3d(*jt);
                    if (lines_intersect3d(A, B)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    // solve direct kinematics using DH parameter table
    virtual std::vector <std::pair <point <real>, point <real>>> dir_kine (
        const point <real> configuration
    ) const override {
        if (configuration.get_dimension() != dh.size()) {
            throw std::length_error("More parameters than links in manipulator");
        }   
        
        for (size_t dim = 0; dim < configuration.get_dimension(); ++dim) {
            const real eps = 1e-6;
            if (configuration[dim] < arm<real>::joint_limits[dim].first - eps || configuration[dim] > arm<real>::joint_limits[dim].second + eps) {
                throw std::domain_error("Configuration coordinates lie outside joint limits");
            }
        } 

        std::vector <real> comp = configuration.get_components();
        std::vector <matrix <real>> homtr(comp.size());
        homtr[0] = dh[0].get_matrix(std::vector<real>{comp[0]});

        // find homogeneus transformation matrix for every segment
        for (size_t lid = 1; lid < comp.size(); ++lid) {
            matrix <real> tmat = dh[lid].get_matrix(std::vector <real>{comp[lid]});
            homtr[lid] = homtr[lid-1] * tmat;
        }

        // multiply each matrix by (0,0,0,1)^T
        // to obtain the coordinates of the link connections
        std::vector <point <real>> connections(comp.size() + 1);
        connections[0] = arm<real>::base;
        matrix <real> zero(point <real>(std::vector <real>{0, 0, 0}));
        
        for (size_t lid = 0; lid < comp.size(); ++lid) {
            matrix <real> origin = homtr[lid] * zero;
            connections[lid + 1] = origin.to_point();
        }

        // connect link contacts to obtain list of links
        std::vector <std::pair <point <real>, point <real>>> links(comp.size());

        for (size_t lid = 0; lid < comp.size(); ++lid) {
            links[lid] = {connections[lid], connections[lid+1]};
        }

        return links;
    }
};


// class implementing antropomorphic arm
template <typename real>
class antropomorphic_arm : public arm3d <real> {
public:
    antropomorphic_arm () {}

    antropomorphic_arm (const antropomorphic_arm <real>& robot) : arm3d <real>(robot.base, robot.link_lengths, robot.joint_limits, robot.dh) {}

    antropomorphic_arm (
        const point <real> _base,
        const std::vector <real>& _link_lengths,
        const std::vector <std::pair <real,real>>& _joint_limits
    ) : arm3d<real>(_base, _link_lengths, _joint_limits) {
        
        if (_link_lengths.size() != 3) {
            throw std::length_error("Antropomorphic arm has 3 links");
        }

        // fill DH parameters for antropomorphic arm, based on link lengths
        // first joint  : {a = 0, d = L1, phi = var, alpha = pi/2}
        // second joint : {a = L2, d = 0, phi = var, alpha = 0}
        // third joint  : {a = L3, d = 0, phi = var, alpha = 0}
        std::vector <std::vector <real>> dhp{
            {0, _link_lengths[0], 0, M_PI_2},
            {_link_lengths[1], 0, 0, 0},
            {_link_lengths[2], 0, 0, 0}
        };

        arm3d <real>::dh.resize(3);
        for (size_t i = 0; i < 3; ++i) {
            arm3d <real>::dh[i] = denavit_hartenberg <real>(dhp[i], "revolute");
        }
    }
};


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

        return false;
    }

    // if robot intersects self
    bool intersects_self (
        const point <real>& configuration
    ) const override {
        auto links = dir_kine(configuration);
    
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
            const real eps = 1e-6;
            if (configuration[dim] < arm<real>::joint_limits[dim].first - eps || configuration[dim] > arm<real>::joint_limits[dim].second + eps) {
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
    std::shared_ptr <arm <real>> robot;

public:
    workspace () {}

    workspace (const workspace <real>& W) : robot(W.robot) {}

    workspace (
        arm <real>* const _robot
    ) : robot(_robot) {}

    workspace (
        const std::shared_ptr <arm <real>>& _robot
    ) : robot(_robot) {}

    // return link lengths
    std::vector <real> get_link_lengths () const {
        return robot->get_link_lengths();
    }

    arm <real>& get_robot () const {
        return *robot;
    }
};


// define 2D workspace
template <typename real>
class workspace2d : public workspace <real> {
private:
    std::vector <polygon<real>> obstacles;

public:
    workspace2d () {}

    workspace2d (const workspace2d <real>& ws) : obstacles(ws.obstacles), workspace <real>(ws.robot) {}

    workspace2d (
        const std::vector <polygon<real>> _obstacles,
        const std::shared_ptr <arm <real>>& _robot
    ) : obstacles(_obstacles), workspace <real>(_robot) {}

    workspace2d (
        const std::vector <polygon<real>> _obstacles,
        arm <real>* const _robot
    ) : obstacles(_obstacles), workspace <real>(_robot) {}

    // test if specific configuration colides
    bool collides (const point <real> configuration) {
        for (const auto& obstacle : obstacles) {
            if (workspace<real>::robot->intersects(obstacle, configuration)) {
                return true;
            }
        }
        if (workspace<real>::robot->intersects_self(configuration)) {
            return true;
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
private:
    std::vector <box <real>> obstacles;

public:
    workspace3d () {}

    workspace3d (const workspace3d <real>& ws) : obstacles(ws.obstacles), workspace <real>(ws.robot) {}

    workspace3d (
        const std::vector <box<real>> _obstacles,
        const std::shared_ptr <arm <real>>& _robot
    ) : obstacles(_obstacles), workspace <real>(_robot) {}

    workspace3d (
        const std::vector <box<real>> _obstacles,
        arm <real>* const _robot
    ) : obstacles(_obstacles), workspace <real>(_robot) {}

    // test if specific configuration colides
    bool collides (const point <real> configuration) {
        for (const auto& obstacle : obstacles) {
            if (workspace<real>::robot->intersects(obstacle, configuration)) {
                return true;
            }
        }
        if (workspace<real>::robot->intersects_self(configuration)) {
            return true;
        }
        return false;
    }
    
    const std::vector <box<real>> get_obstacles () const {
        return obstacles;
    }
};