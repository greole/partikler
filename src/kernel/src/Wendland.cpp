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

#include "Wendland.hpp"

#include "ObjectRegistry.hpp"

// TODO Implement 2d and 3d kernel by deriving from base class

Wendland::Wendland(
    const std::string &model_name,
    YAML::Node parameter,
    ObjectRegistry &objReg,
    float hfact)
    : Model(model_name, parameter, objReg),
      h_(read_or_default_coeff<Scalar>("h", 1.0)), ih_(1.0 / h_),
      W_fak2_(hfact / (h_ * h_)), dW_fak2_(hfact / (h_ * h_ * h_)),
      pos_(objReg.get_pos()),
      np_(objReg.get_object<Field<std::vector<NeighbourPair>>>(
          "neighbour_pairs")),
      W_(objReg.create_field<ScalarField>("KernelW")),
      dWdx_(objReg.create_field<KernelGradientField>("KerneldWdx")) {
    objReg.create_generic<Scalar>("h", h_);
}

void Wendland::execute() {

    // log().set_scope("Wendland::execute()");
    log().info_begin() << "Computing Kernel";

    const size_t size {np_.size()};

    // Resize kernel
    W_.resize(size);
    dWdx_.resize(size);

    for (size_t pid = 0; pid < size; pid++) {

        auto [oid, nid] = np_[pid];

        Vec3 const &opos = pos_[oid];
        Vec3 const &npos = pos_[nid];

        // auto [len, lenVo, lenVn] = sd_[pid];

        const auto len = opos - npos;
        Scalar const maglen = std::sqrt(squared_length(len));

        const Scalar q {maglen * ih_};

        if (q > 2.) {
            log().warn() << "Outside kernel radius";
            W_[pid] = 0.0;
            dWdx_[pid] = {0, 0, 0};
            continue;
        }

        const Scalar q3 = (q - 2.);
        const Scalar qfac2 = q3 * q3;
        const Scalar qfac4 = qfac2 * qfac2;

        Scalar q2 = 2. * q;
        q2 += 1.;

        W_[pid] = qfac4 * q2 * W_fak2_;

        const Scalar prefact = 10. * qfac2 * q * dW_fak2_;
        if (maglen > 0.0) {
            for (int j = 0; j < 3; j++) {
                dWdx_[pid][j] = (Scalar)len[j] / maglen * prefact;
            }
        } else {
            for (int j = 0; j < 3; j++) {
                dWdx_[pid] = {0.0, 0.0, 0.0};
            }
            log().warn() << "Neighbour sum == 0";
        }
    }

    log().info_end();
}

Wendland2D::Wendland2D(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Wendland(model_name, parameter, objReg, 7. / 64.) {}

Wendland3D::Wendland3D(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Wendland(
          model_name,
          parameter,
          objReg,
          21. / (256. *
                 read_or_default_coeff_impl<Scalar>(parameter, "h", 1.0))) {}

REGISTER_DEF_TYPE(KERNEL, Wendland2D);
REGISTER_DEF_TYPE(KERNEL, Wendland3D);
