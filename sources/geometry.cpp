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
