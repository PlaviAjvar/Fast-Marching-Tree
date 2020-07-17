#pragma once

#include <vector>
#include <ostream>

// class implementing point
template <typename real>
class point {
private:
    std::vector <real> components;

public:
    point () {}

    point (const size_t dim) {
        (this->components).resize(dim);
    }

    point (const std::vector <real>& comp) : components(comp) {}

    point (const point <real>& P) : components(P.get_components()) {}

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
};

// operator overload for outputing point<double>

std::ostream& operator << (std::ostream& os, const point <double> P);

// euclidean distance
template <typename real>
real euclidean_distance (const point<real>& A, const point<real>& B) {
    if (A.get_dimension() != B.get_dimension()) {
        throw std::length_error("Dimensions of points differ");
    }

    size_t dimension = A.get_dimension();
    real total = 0;

    for (int index = 0; index < dimension; ++index) {
        total += (A[index] - B[index]) * (A[index] - B[index]);
    }

    return total;
}