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

#ifndef Solenthaler_H
#define Solenthaler_H

#include <string> // for string

#include "Equation.hpp"
#include "Field.hpp" // for ScalarField, VectorField
#include "FieldOps.hpp"
#include "Models.hpp" // for ScalarFieldEquation, ModelRegister (ptr only)
#include "Scalar.hpp"
#include "SearchCubes.hpp"
#include "yaml-cpp/yaml.h"

#include <boost/yap/print.hpp>
#include <boost/yap/yap.hpp>

#include "Scalar.hpp"
class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class SolenthalerGradient : public ScalarGradientEquation {

    REGISTER_DEC_TYPE(SolenthalerGradient);

  private:
    // In
    // Density
    ScalarFieldEquation &conti_;
    ScalarFieldEquation &pressure_;

    // Coeffs
    const Scalar rho_0_;
    const Scalar eta_;
    const int min_iter_;
    const ScalarField &mp_;

  public:
    SolenthalerGradient(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
