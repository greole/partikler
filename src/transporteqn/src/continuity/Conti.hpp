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

#ifndef CONTI_H
#define CONTI_H

#include <string> // for string

#include "Equation.hpp"
#include "Field.hpp" // for ScalarField, PointField
#include "FieldOps.hpp"
#include "Models.hpp" // for ScalarFieldEquation, ModelRegister (ptr only)
#include "Scalar.hpp"
#include "SearchCubes.hpp"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class ParticleMass : public ScalarFieldEquation {

    REGISTER_DEC_TYPE(ParticleMass);

  private:
    // Coeffs

  public:
    ParticleMass(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

class Conti : public ScalarFieldEquation {

    REGISTER_DEC_TYPE(Conti);

  private:
    // Coeffs
    const Scalar rho_0_;

    ScalarField &mp_;

  public:
    Conti(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

class TransientConti : public ScalarFieldEquation {

    REGISTER_DEC_TYPE(TransientConti);

  private:
    // Coeffs

    VectorField &u_;

    ScalarField &mp_;

  public:
    TransientConti(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
