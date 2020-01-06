#include "Field.hpp"

PointField &operator+=(PointField &a, VectorField &b) {
    for (size_t i = 0; i < a.size(); i++) {
        a[i] += {b[i][0], b[i][1], b[i][2]};
    }
    return a;
}
