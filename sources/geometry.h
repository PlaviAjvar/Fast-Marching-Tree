#pragma once

#include <vector>
#include <ostream>
#include <functional>
#include <iostream>

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

    real getx () {
        return point<real>::components[0];
    }

    real gety () {
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

    real getz() const {
        return point<real>::components[2];
    }
};


// arbitrary polygon in 2d space
template <typename real>
class polygon {
private:
    // given in clockwise order
    std::vector <std::pair <point2d <real>, point2d <real>>> edges;

public:
    polygon () {}

    polygon (const polygon& P) : edges(P.edges) {}

    polygon (const std::vector <point2d <real>>& _edges) : edges(_edges) {}

    bool in_polygon (const point2d <real> point);
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

    bool in_box (const point <real> point) const {
        if (box_limits.size() != point.get_dimension()) {
            throw std::length_error("Dimensionality of box differs from dimensionality of point");
        }

        size_t dimension = point.get_dimension();

        for (size_t dim = 0; dim < dimension; ++dim) {
            // coordinate falls outside the feasible range
            if (point[dim] < box_limits[dim].first || point[dim] > box_limits[dim].second) {
                return true;
            }
        }

        return false;
    }

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
                // for (real z = box_limits[2].first; z <= box_limits[2].second; z += stepsize) {
                //     xr.push_back(x);
                //     yr.push_back(y);
                //     zr.push_back(z);
                // }

                xr.push_back(x);
                yr.push_back(y);
                zr.push_back(0);
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

    for (size_t index = dimension-2; index >= 0; --index) {
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
