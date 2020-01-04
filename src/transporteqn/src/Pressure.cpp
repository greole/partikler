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

# include "Pressure.hpp"

Pressure::Pressure(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : FloatFieldEquation(model_name, parameter, objReg, objReg.create_field<FloatField>("p", p_0_)),
      c_(read_or_default_coeff<float>("c", 300.0)),
      rho_0_(read_or_default_coeff<float>("rho_0", 1.0)),
      gamma_(read_or_default_coeff<float>("gamma", 1.4)),
      p_0_(read_or_default_coeff<float>("p_0", 10000)),
      prefac_(c_ * c_ * rho_0_ / gamma_),
      rho_(objReg.get_object<FloatField>("rho")),
      p_(objReg.create_field<FloatField>("p", p_0_)),
      dp_(objReg.create_field<VectorField>(
          "dp", zero<VectorField::value_type>::val, {"dpx", "dpy", "dpz"})) {};

void Pressure::execute() {

    log().info_begin() << "Computing pressure";

    // const ScalarField tmp0 = (rho_/rho_0_).pow(gamma_);
    // const ScalarField tmp1 = ((tmp0 - 1.0)*prefac_)+p_0_;


    // assign(p,
    //     prefac_ * ((rho_ / rho_0_).pow(gamma_) - 1.0) + p_0_,
    //     p_
    // );

    // p_ = tmp1;
    // // log().info_end();

    // // log().info_begin() << "Computing gradient";

    // const FloatField prho = p_/(rho_*rho_);
    // const FloatField tmp_ab = prho.add_ab(np_);

    // // reset
    // dp_.set_uniform({0,0,0});

    // dp_.weighted_sum(np_, tmp_ab, dW_);

    log().info_end();
};

REGISTER_DEF_TYPE(TRANSPORTEQN, Pressure);
