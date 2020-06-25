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

#ifndef PARTIKLER_GRAVITY_INCLUDED
#define PARTIKLER_GRAVITY_INCLUDED

#include <string> // for string

#include "Equation.hpp"
#include "Field.hpp" // for ScalarField, PointField
#include "FieldOps.hpp"
#include "Models.hpp" // for ScalarFieldEquation, ModelRegister (ptr only)
#include "SearchCubes.hpp"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class Gravity : public VectorFieldEquation {

    REGISTER_DEC_TYPE(Gravity);

  private:
    // Coeffs
    const Vec3 gravity_;

  public:
    Gravity(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
