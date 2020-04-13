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

#include "Gauss.hpp"

#include "ObjectRegistry.hpp"
#include "math.h"

// TODO Implement 2d and 3d kernel by deriving from base class

Gauss::Gauss(
    const std::string &model_name,
    YAML::Node parameter,
    ObjectRegistry &objReg,
    Scalar hfact)
    : Model("KernelEqn", parameter, objReg),
      h_(read_or_default_coeff<Scalar>("h", 1.0)), ih_(1.0 / h_),
      alpha_(hfact / (h_ * h_)), pos_(objReg.get_pos()),
      np_(objReg.get_object<Field<std::vector<NeighbourPair>>>(
          "neighbour_pairs")),
      W_(objReg.create_field<ScalarField>("KernelW")),
      dWdx_(objReg.create_field<KernelGradientField>("KerneldWdx")) {
    objReg.create_generic<Scalar>("h", h_);
    objReg.reference_clone("KerneldWdx", "KerneldWdxNeighbour");
}

void Gauss::execute() {

    // log().set_scope("Gauss::execute()");
    log().info_begin() << "Computing Kernel";

    const size_t size {np_.size()};

    // Resize kernel
    W_.resize(size);
    dWdx_.resize(size);

    for (size_t pid = 0; pid < size; pid++) {

        auto [oid, nid] = np_[pid];

        Vec3 const &opos = pos_[oid];
        Vec3 const &npos = pos_[nid];

        const auto len = opos - npos;
        Scalar const maglen = std::sqrt(squared_length(len));

        const Scalar q {maglen * ih_};

        if (q > 2.) {
            log().warn() << "Outside kernel radius";
            W_[pid] = 0.0;
            dWdx_[pid] = {0, 0, 0};
            continue;
        }

        Scalar W = exp(-q * q) * alpha_;

        W_[pid] = W;

        const Scalar prefact = -2.0 * q * W;

        if (maglen > 0.0) {
            for (int j = 0; j < 3; j++) {
                dWdx_[pid][j] = (Scalar)len[j] / maglen * prefact;
            }
        } else {
            dWdx_[pid] = {0.0, 0.0, 0.0};
            log().warn() << "Neighbour sum == 0";
        }
    }

    log().info_end();
}

Gauss2D::Gauss2D(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Gauss(model_name, parameter, objReg, 1.0 / M_PI) {}

Gauss3D::Gauss3D(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Gauss(
          model_name,
          parameter,
          objReg,
          1.0 / ((pow(M_PI, 3 / 2)) *
                 read_or_default_coeff_impl<Scalar>(parameter, "h", 1.0))) {}

REGISTER_DEF_TYPE(KERNEL, Gauss2D);
REGISTER_DEF_TYPE(KERNEL, Gauss3D);
