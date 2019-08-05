/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

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
