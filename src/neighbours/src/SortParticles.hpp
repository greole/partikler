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

#ifndef SORTPARTRICLES_H
#define SORTPARTRICLES_H

#include <string> // for string
#include <vector> // for vector

#include "Field.hpp"  // for Field (ptr only), PointField, SizeTField
#include "Models.hpp" // for Model, ModelRegister (ptr only), REGISTER_DEC_...
#include "SearchCubes.hpp"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML
struct SearchCube;
struct SearchCubeDomain;
template <class T> class Generic;

class CountingSortParticles : public Model {

    REGISTER_DEC_TYPE(CountingSortParticles);

  private:
    // In
    PointField &pos_;

    Field<std::vector<SearchCube>> &sc_;

    //  Sorting indexes
    SizeTField &si_;

    Generic<SearchCubeDomain> &scd_;

  public:
    CountingSortParticles(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();

    void reorder_fields();
};

class CountingSortParticlesVec3 : public Model {

    REGISTER_DEC_TYPE(CountingSortParticlesVec3);

  private:
    // In
    VectorField &pos_;

    Field<std::vector<SearchCube>> &sc_;

    //  Sorting indexes
    SizeTField &si_;

    Generic<SearchCubeDomain> &scd_;

  public:
    CountingSortParticlesVec3(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();

    void reorder_fields();
};

#endif
