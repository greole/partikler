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
#include "Conti.hpp"

#include "Time.hpp"

Momentum::Momentum(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : VectorFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "u", zero<VectorField::value_type>::val, {"U", "V", "W"})),
      conti_(objReg.get_or_create_model<Conti>("Conti", parameter, objReg)),
      tau_(objReg.get_object<VectorFieldEquation>("Viscosity")),
      p_(objReg.get_object<ScalarFieldEquation>("Pressure")),
      g_(objReg.get_object<VectorFieldEquation>("Gravity")),
      fc_(objReg.get_object<VectorFieldEquation>("CohesionAkinci"))
{}

void Momentum::execute() {

    log().info_begin() << "Computing du/dt";

    int it = time_.get_current_timestep();

    ScalarField &rho = conti_.get(time_.get_current_timestep());

    solve(df_, tau_.get(it) / rho - p_.get_dx(it) / rho + g_.get(it) + fc_.get(it));

    log().info_end();

    log().info_begin() << "Computing velocity";


    for (size_t i = 0; i < f_.size(); i++) {
        // f_[i] = f_[i] / 1.3;
        f_[i] += df_[i] * time_.get_deltaT();
    }

    log().info_end();

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Momentum);
