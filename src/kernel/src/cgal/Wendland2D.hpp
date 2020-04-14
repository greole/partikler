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

#include "Scalar.hpp"

#include "Field.hpp"       // for Field (ptr only), ScalarField, KernelGradi...
#include "Models.hpp"      // for Model, ModelRegister (ptr only), REGISTER...
#include "SearchCubes.hpp" // for NeighbourFieldAB
#include "cgal/CGALHelper.hpp"
#include "cgal/CGALField.hpp"
#include "yaml-cpp/yaml.h"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML
struct STLSurfaceDist;

class STLWendland2D : public Model {

    REGISTER_DEC_TYPE(STLWendland2D);

  private:
    // Coeffs
    const Scalar h_;       // Smoothing length
    const Scalar ih_;      // Inverse Smoothing length
    const Scalar W_fak2_;  // = 7. / (64. * M_PI * h * h);
    const Scalar dW_fak2_; // = 7. / (64. * M_PI * h * h * h);

    const NeighbourFieldAB &np_;
    const Field<std::vector<STLSurfaceDist>> &sd_;

    // Out
    // Kernel &kernel                               // Kernel field
    ScalarField &W_;

    KernelGradientField &dWdx_;

    // KernelGradient as seen from the neighbouring particle
    // since KernelGradients on the surfaces are not isotropic
    KernelGradientField &dWdxn_;

  public:
    STLWendland2D(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
