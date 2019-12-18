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

#include "Models.hpp"
#include "SearchCubes.hpp"
#include "cgal/CGALHelper.hpp"

#include "yaml-cpp/yaml.h"

class STLWendland2D : public Model {

    REGISTER_DEC_TYPE(STLWendland2D);

  private:
    // Coeffs
    const float h_;  // Smoothing length
    const float ih_; // Inverse Smoothing length
    const float W_fak2_;  // = 7. / (64. * M_PI * h * h);
    const float dW_fak2_; // = 7. / (64. * M_PI * h * h * h);

    // In
    const PointField &pos_; // Particle positions

    const NeighbourFieldAB &np_;
    const Field<std::vector<STLSurfaceDist>> &sd_;

    // Out
    // Kernel &kernel                               // Kernel field
    FloatField &W_;
    KernelGradientField &dWdx_;

  public:
    STLWendland2D(
        const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg);

    void execute();
};

#endif
