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

#include "Cole.hpp"

#include "Time.hpp"

Cole::Cole(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : ScalarFieldEquation(
          "Pressure",
          parameter,
          objReg,
          objReg.create_field<ScalarField>("p", 10000)),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      c_(read_or_default_coeff<Scalar>("c", 300.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      gamma_(read_or_default_coeff<Scalar>("gamma", 1.4)),
      p_0_(read_or_default_coeff<Scalar>("p_0", 10000)),
      prefac_(c_ * c_ * rho_0_ / gamma_) {
    maxDt_ = 0.5 / c_ * h_;
    time_.set_model_timestep(model_name, maxDt_);
    init_limits();
}

void Cole::execute() {

    auto &rho = conti_.get();

    auto pow_s = Pow_Wrapper<Scalar>(gamma_);
    auto pow = boost::yap::make_terminal(pow_s);

    solve(prefac_ * (pow(rho / rho_0_) - 1.0) + p_0_);

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(PRESSURE, Cole);
