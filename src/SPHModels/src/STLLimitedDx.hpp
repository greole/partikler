#ifndef STLLIMITEDDX_H
#define STLLIMITEDDX_H

#include "SPHDatastructures.hpp"

SPHVectorField STL_limited_dx(
    SPHVectorField &u,
    float dt,
    // const std::vector<Point> &opoints,
    SPHField<Facet_handle> &facets,
    const SPHIntField &type,
    const SPHSizeTField &idx,
    const SPHPointField &pos);

#endif
