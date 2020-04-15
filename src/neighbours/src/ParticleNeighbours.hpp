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

#ifndef PARTIKLER_PARTICLENEIGHBOURS_INCLUDED
#define PARTIKLER_PARTICLENEIGHBOURS_INCLUDED

#include <string> // for string
#include <vector> // for vector

#include "Field.hpp"             // for FieldAB, Field (ptr only), PointField
#include "Models.hpp"            // for Model, ModelRegister (ptr only)
#include "SearchCubes.hpp"       // for SearchCube, NeighbourFieldAB, Searc...

#include "Vec3.hpp"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML
template <class T> class Generic;

class SPHParticleNeighbours : public Model {

    REGISTER_DEC_TYPE(SPHParticleNeighbours);

    using SearchCubeFieldAB = FieldAB<Field<std::vector<SearchCube>>>;
    using DistAB = FieldAB<Field<std::vector<Vec3>>>;

  private:
    // Coeffs
    float dx_;

    // In
    VectorField &pos_;

    SearchCubeFieldAB &sc_;

    // Out
    NeighbourFieldAB &np_;

    DistAB &sd_;

    // Regular data member
    Generic<SearchCubeDomain> &scd_;

    float search_cube_size_ = 1.0;

  public:
    SPHParticleNeighbours(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();

    void update_search_cube_domain();
};

#endif
