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

#include "Kernel.hpp"

STLWendland2D::STLWendland2D(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)
    : SPHModel(model_name, parameter, runTime),
      h_(read_or_default_coeff<float>("h", 1.0)), ih_(1.0 / h_),
      W_fak2_(7. / (64. * M_PI * h_ * h_)),
      dW_fak2_(7. / (64. * M_PI * h_ * h_ * h_)),
      pos_(get_runTime().get_particle_positions()),
      np_(get_runTime().get_obj<SPHField<searchcubes::NeighbourPair>>(
          "neighbour_pairs")),
      sd_(get_runTime().get_obj<SPHField<STLSurfaceDist>>("surface_dist")),
      W_(get_runTime().create_field<SPHFloatField>("KernelW")),
      dWdx_(get_runTime().create_field<SPHField<VectorPair>>("KerneldWdx")) {};

void STLWendland2D::execute() {

    // log().set_scope("STLWendland2D::execute()");
    log().info_begin() << "Computing Kernel";

    const size_t size {np_.size()};

    // Resize kernel
    W_.resize(size);
    dWdx_.resize(size);

    for (size_t pid = 0; pid < size; pid++) {

        auto [oid, nid] = np_[pid];

        const Point &opos = pos_[oid];
        const Point &npos = pos_[nid];

        auto [len, lenVo, lenVn] = sd_[pid];

        const float q {len * ih_};

        if (q > 2.) {
            log().warn() << "Outside kernel radius";
            W_[pid] = 0.0;
            dWdx_[pid] = {{0,0,0}, {0,0,0}};
            continue;
        }

        const float q3 = (q - 2.);
        const float qfac2 = q3 * q3;
        const float qfac4 = qfac2 * qfac2;

        float q2 = 2. * q;
        q2 += 1.;

        W_[pid] = qfac4 * q2 * W_fak2_;

        const float prefact = 10. * qfac2 * q * dW_fak2_;
        if( len != 0.0) {
          for(int j=0;j<3;j++) {
              dWdx_[pid].on[j] = (float)lenVo[j] / len * prefact;
              dWdx_[pid].no[j] = (float)lenVn[j] / len * prefact;
          }
        } else {
          for(int j=0;j<3;j++) {
              dWdx_[pid] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
          }
          log().warn() << "Neighbour sum == 0";
        }
    }

    log().info_end();
};

REGISTER_DEF_TYPE(KERNEL, STLWendland2D);
