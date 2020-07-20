#include "geometry.h"
 
std::ostream& operator << (std::ostream& os, const point <double> P) {
    std::vector <double> comp = P.get_components();
    os << "(";

    for (std::vector<double>::iterator it = comp.begin(); it != comp.end(); ++it) {
        if (it != comp.begin()) os << ", ";
        os << *it;
    }

    os << ")";
    return os;
}

// walk along direction from nearest to sample

point <double> walk (
    const point <double>& nearest, 
    const point <double>& sample, 
    const double scaler
) {

    point <double> diff(sample - nearest);
    return nearest + diff * scaler;
}

// check if path is clear between two points

bool path_clear(
    const point <double> A,
    const point <double> B,
    const std::function <bool(point<double>)>& collision_check,
    const double stepsize,
    const double epsilon
) {

    for (double scaler = stepsize; scaler <= 1 + epsilon; scaler += stepsize) {
        // point 1 step further down line
        point <double> step_further = walk(A, B, scaler);

        // if no collision update new configuration
        if (collision_check(step_further)) {
            return false;
        }
    }
    
    return true;
}