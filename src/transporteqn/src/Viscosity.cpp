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
#include "Viscosity.hpp"

Viscosity::Viscosity(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "tau",
              zero<VectorField::value_type>::val,
              {"taux", "tauy", "tauz"})
              ),
      conti_(objReg.get_or_create_model<Conti>("Conti", parameter, objReg)),
      nu_(read_or_default_coeff<float>("nu", 1e-05)),
      u_(objReg.create_field<VectorField>(
          "u", zero<VectorField::value_type>::val, {"U", "V", "W"})),
      pos_(objReg.get_object<PointField>("Pos"))
{}

void Viscosity::execute() {

    log().info_begin() << "Computing dnu";

    // auto& u = momentum_.get();
    FloatField& rho = conti_.get(time_.get_current_timestep());

    // sum_AB_dW(f_, np_, dW_,
    //       ( A(u_) + B(u_) )/( A(rho) + B(rho) )
    //           *  ((ab(u_) * ab(pos_))/ (ab(pos_))*(ab(pos_))
    //     );


    // const SPHFloatField tmp0 = nu.add_ab(pn) / rho.mult_ab(pn);

    // const VectorField dxp = particle_distance_vec(pos_, np_);

    // // const SPHFloatField tmp1 = (u.sub_ab(pn) * dxp) / dxp.norm();

    // // const SPHFloatField tmp = tmp0*tmp1;

    // const FloatField tmp = (u_.sub_ab(np_) * dxp)/dxp.norm();

    // // TODO Reset
    // dnu_.set(Vector {0,0,0});

    // // weighted sum
    // const size_t size = pn.ids.size();
    // for (size_t i = 0; i < size; i++) {
    //     size_t ownId = pn.ids[i].ownId;
    //     size_t neighId = pn.ids[i].neighId;

    //     for (int j = 0; j < 3; j++) {
    //         dnu[ownId][j]   -= tmp[i]  * kernel.dWdxo[i][j];
    //         dnu[neighId][j] -= tmp[i]  * kernel.dWdxn[i][j];
    //     }
    // }

    // dnu_.weighted_sum(np_, tmp, dW_);

    // dnu_ *= nu_;

    log().info_end();
    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Viscosity);
