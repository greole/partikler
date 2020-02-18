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

#include "Momentum.hpp"
#include "Equation.hpp"

#include "Time.hpp"

Momentum::Momentum(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : VectorFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.velocity()),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      dp_(objReg.get_object<ScalarGradientEquation>("PressureGradient")),
      dtau_(objReg.get_object<VectorGradientEquation>("Viscosity")),
      p_(objReg.get_object<ScalarFieldEquation>("Pressure")),
      g_(objReg.get_object<VectorFieldEquation>("Gravity")),
      fc_(objReg.get_object<VectorFieldEquation>("Akinci")) {}

void Momentum::execute() {

    log().info_begin() << "Computing du/dt";

    int it = time_.get_current_timestep();

    auto &rho = conti_.get(it);
    auto &dp = dp_.get(it);
    auto &dtau = dtau_.get(it);

    size_t a2 = 0;

    auto ddt = boost::yap::make_terminal(ddt_);

    solve(ddt(dtau / rho - dp / rho + g_.get(it) + fc_.get(it)));

    log().info_end();

    Scalar maxAbsU = 0;
    for(auto &el: f_) {
        Scalar absU = mag(el);
        if (absU > maxAbsU) maxAbsU = absU;
    }
    if (maxAbsU > 0) {
        maxDt_ = 0.5*h_/maxAbsU;
        time_.set_model_timestep("Momentum", maxDt_);
    }

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Momentum);
