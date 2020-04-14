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

#include <vector> // for vector

#include "Vec3.hpp"                         // for Vec3
#include "cgal/CGALHelper.hpp"              // for inwardPointingEdgeNormal
#include <CGAL/Cartesian/Vector_3.h>        // for operator==
#include <CGAL/HalfedgeDS_iterator.h>       // for I_HalfedgeDS_facet_circ
#include <CGAL/HalfedgeDS_list.h>           // for HalfedgeDS_in_place_list...
#include <CGAL/In_place_list.h>             // for In_place_list_iterator
#include <CGAL/Kernel/global_functions_3.h> // for operator*, operator==
#include <CGAL/Line_3.h>                    // for Line_3<>::Point_3
#include <CGAL/Point_3.h>                   // for Point_3
#include <CGAL/Simple_cartesian.h>          // for Simple_cartesian, Cartes...
#include <CGAL/Vector_3.h>                  // for Vector_3
#include <array>                            // for array
#include <ext/alloc_traits.h>               // for __alloc_traits<>::value_...
#include <math.h>                           // for sqrt
#include <stddef.h>                         // for size_t

#include "Field.hpp" // for VectorField, Field (ptr only), IntF...
#include "cgal/CGALHelper.hpp"
#include "cgal/CGALTYPEDEFS.hpp" // for Facet_handle
#include "cgal/CGALField.hpp" // for Facet_handle

void STL_limited_dx(
    VectorField &u,
    float dt,
    Field<std::vector<Facet_handle>> &facets,
    const IntField &type,
    const SizeTField &idx,
    const PointField &pos,
    VectorField &ret);

#endif
