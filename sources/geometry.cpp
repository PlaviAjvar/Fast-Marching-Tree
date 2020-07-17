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