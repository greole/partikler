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

#ifndef PARTIKLER_FIXEDVALUE_INCLUDED_H
#define PARTIKLER_FIXEDVALUE_INCLUDED_H

#include <map>    // for map
#include <string> // for string

#include "Models.hpp" // for Model, ModelRegister (ptr only), REGISTER_DEC_...
#include "Scalar.hpp" // for Model, ModelRegister (ptr only), REGISTER_DEC_...
#include "Vec3.hpp"   // for Model, ModelRegister (ptr only), REGISTER_DEC_...
#include "yaml-cpp/yaml.h"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class FixedValue : public Model {

    REGISTER_DEC_TYPE(FixedValue);

  private:
    // maps field and boundary name to fixed uniform values
    std::map<std::string, std::map<std::string, Scalar>> float_fields_ {};
    std::map<std::string, std::map<std::string, Vec3>> vec_fields_ {};

  public:
    FixedValue(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
