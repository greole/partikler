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

#ifndef STLPOSINTEGRATOR_H
#define STLPOSINTEGRATOR_H

#include <string> // for string
#include <vector> // for vector

#include "Field.hpp" // for Field (ptr only), IntField, PointField
#include "FieldOps.hpp"
#include "Models.hpp"            // for Model, ModelRegister (ptr only)
#include "cgal/CGALTYPEDEFS.hpp" // for Facet_handle
#include "stl/STLLimitedDx.hpp"
#include "yaml-cpp/yaml.h"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class STLPosIntegrator : public Model {

    REGISTER_DEC_TYPE(STLPosIntegrator);

  private:
    // const std::vector<Point> opoints_;
    const IntField &type_;
    Field<std::vector<Facet_handle>> &facets_;
    const SizeTField &idx_;

    // Out
    VectorField &u_;
    PointField &pos_;
    TimeGraph &time_;

  public:
    STLPosIntegrator(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
