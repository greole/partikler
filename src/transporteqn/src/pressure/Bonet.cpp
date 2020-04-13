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

#include "Bonet.hpp"

#include "Time.hpp"

BonetGradient::BonetGradient(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : ScalarGradientEquation(
          "PressureGradient",
          parameter,
          objReg,
          objReg.get_object<ScalarField>("p")),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      pressure_(objReg.get_object<ScalarFieldEquation>("Pressure")),
      mp_(objReg.get_object<ScalarField>("mp")) {}

void BonetGradient::execute() {

    auto &p = pressure_.get();
    auto &rho = conti_.get();

    auto &dW = this->get_objReg().template get_object<KernelGradientField>(
        "KerneldWdx");

    auto sum_AB_e = Sum_AB_dW_sym<VectorField>(f_, np_, dW, dW);
    auto sum_AB_dW_e = boost::yap::make_terminal(sum_AB_e);

    solve(1.0 / rho * sum_AB_dW_e(mp_.b() / rho.b() * (p.a() + p.b())));

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(PRESSURE, BonetGradient);
