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

#include "Pressure.hpp"

Pressure::Pressure(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : FloatFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<FloatField>("p", 10000)),
      conti_(objReg.get_or_create_model<Conti>("Conti", parameter, objReg)),
      c_(read_or_default_coeff<float>("c", 300.0)),
      rho_0_(read_or_default_coeff<float>("rho_0", 1.0)),
      gamma_(read_or_default_coeff<float>("gamma", 1.4)),
      p_0_(read_or_default_coeff<float>("p_0", 10000)),
      prefac_(c_ * c_ * rho_0_ / gamma_),
      p(f_), dp(df_)
{}

void Pressure::execute() {


    auto& rho = conti_.get(time_.get_current_timestep());

    auto pow = boost::yap::make_terminal(Pow_Wrapper<float>(gamma_));

    // auto pow = Pow<float>(gamma_);

    log().info_begin() << "Computing pressure";

    solve(p, prefac_ * (  pow(rho/rho_0_ ) - 1.0) + p_0_);

    log().info_end();

    log().info_begin() << "Computing gradient";

    sum_AB_dW(dp, np_, dW_, (A(p) + B(p))/B(rho));

    log().info_end();

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Pressure);
