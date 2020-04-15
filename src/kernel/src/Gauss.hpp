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

#ifndef KERNEL_H
#define KERNEL_H

#include <string> // for string
#include <vector> // for vector

#include "Field.hpp"  // for Field (ptr only), ScalarField, PointField
#include "Models.hpp" // for ModelRegister (ptr only), REGISTER_DEC_TYPE
#include "SearchCubes.hpp"
#include "yaml-cpp/yaml.h"

#include "Scalar.hpp"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML
struct NeighbourPair;
struct Vec3;

class Gauss : public Model {

  private:
    // Coeffs
    const Scalar h_;  // Smoothing length
    const Scalar ih_; // Inverse Smoothing length
    // // 3d
    // const float W_fak2 = 21. / (256. * M_PI * h * h * h);
    // const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);
    // 2d
    const Scalar alpha_; // = 7. / (64. * M_PI * h * h);

    // In
    const VectorField &pos_; // Particle positions

    const Field<std::vector<NeighbourPair>> &np_;
    // const Field<STLSurfaceDist> &sd_;

    // Out
    // Kernel &kernel                               // Kernel field
    ScalarField &W_;
    Field<std::vector<Vec3>> &dWdx_;

  public:
    Gauss(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        Scalar hfact);

    void execute();
};

class Gauss2D : public Gauss {

    REGISTER_DEC_TYPE(Gauss2D);

  public:
    Gauss2D(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);
};

class Gauss3D : public Gauss {

    REGISTER_DEC_TYPE(Gauss3D);

  public:
    Gauss3D(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);
};

#endif
