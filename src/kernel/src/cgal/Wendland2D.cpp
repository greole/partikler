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

#include "Wendland2D.hpp"

#include "ObjectRegistry.hpp"

STLWendland2D::STLWendland2D(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      h_(read_or_default_coeff<Scalar>("h", 1.0)), ih_(1.0 / h_),
      W_fak2_(7. / (64. * M_PI * h_ * h_)),
      dW_fak2_(7. / (64. * M_PI * h_ * h_ * h_)),
      np_(objReg.get_object<NeighbourFieldAB>("neighbour_pairs")),
      sd_(objReg.get_object<Field<std::vector<STLSurfaceDist>>>(
          "surface_dist")),
      W_(objReg.create_field<ScalarField>("KernelW")),
      dWdx_(objReg.create_field<KernelGradientField>("KerneldWdx")),
      dWdxn_(objReg.create_field<KernelGradientField>("KerneldWdxNeighbour")) {}

void STLWendland2D::execute() {

    // log().set_scope("STLWendland2D::execute()");
    log().info_begin() << "Computing Kernel";

    const size_t size {np_.size()};

    // Resize kernel
    W_.resize(size);
    dWdx_.resize(size);

    for (size_t pid = 0; pid < size; pid++) {

        auto [len, lenVo, lenVn] = sd_[pid];

        const Scalar q {len * ih_};

        if (q > 2.) {
            log().warn() << "Outside kernel radius";
            W_[pid] = 0.0;
            dWdx_[pid] = {{0, 0, 0}};
            dWdxn_[pid] = {{0, 0, 0}};
            continue;
        }

        const Scalar q3 = (q - 2.);
        const Scalar qfac2 = q3 * q3;
        const Scalar qfac4 = qfac2 * qfac2;

        Scalar q2 = 2. * q;
        q2 += 1.;

        W_[pid] = qfac4 * q2 * W_fak2_;

        const float prefact = 10. * qfac2 * q * dW_fak2_;
        if (len != 0.0) {
            Vec3 lenVoV = {
                (Scalar)lenVo.x(), (Scalar)lenVo.y(), (Scalar)lenVo.z()};
            Vec3 lenVnV = {
                (Scalar)lenVn.x(), (Scalar)lenVn.y(), (Scalar)lenVn.z()};
            dWdx_[pid] = lenVoV / len * prefact;
            dWdxn_[pid] = lenVnV / len * prefact;
        } else {
            dWdx_[pid] = {{0.0, 0.0, 0.0}};
            dWdxn_[pid] = {{0.0, 0.0, 0.0}};
            log().warn() << "Neighbour sum == 0";
        }
    }

    log().info_end();
}

REGISTER_DEF_TYPE(KERNEL, STLWendland2D);
