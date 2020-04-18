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

#include "Conti.hpp"

#include "Time.hpp"
#include <math.h>

ParticleMass::ParticleMass(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ScalarFieldEquation(
          "Conti",
          parameter,
          objReg,
          objReg.create_field<ScalarField>("mp", 0.0)) {

    for (size_t i = 0; i < f_.size(); i++) {
        // get id
        auto &fieldIdMap(objReg.get_object<FieldIdMap>("FieldIdMap"));
        Material m = fieldIdMap.getMaterial(this->id_[i]);
        f_[i] = m.getMp();
    }
}

void ParticleMass::execute() {}

Conti::Conti(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ScalarFieldEquation(
          "Conti",
          parameter,
          objReg,
          objReg.create_field<ScalarField>("rho", 0.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      mp_(objReg.get_object<ScalarField>("mp")) {

    init_limits();
}

void Conti::execute() {

    log().info_begin() << "Computing density";

    auto sab = sum_AB();
    auto sum_AB_o = boost::yap::make_terminal(sab);

    solve(mp_ * sum_AB_o(W_));

    // set iteration
    iteration_ = time_.get_current_timestep();

    log().info_end();
}

TransientConti::TransientConti(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ScalarFieldEquation(
          "Conti",
          parameter,
          objReg,
          objReg.create_field<ScalarField>(
              "rho",
              read_or_default_coeff_impl<Scalar>(parameter, "rho_0", 1.0))),
      u_(objReg.velocity()), mp_(objReg.get_object<ScalarField>("mp")) {

    // init field from material properties
    for (size_t i = 0; i < f_.size(); i++) {
        // get id
        auto &fieldIdMap(objReg.get_object<FieldIdMap>("FieldIdMap"));
        Material m = fieldIdMap.getMaterial(this->id_[i]);
        f_[i] = m.getRho();
    }

    init_limits();
}

void TransientConti::execute() {

    auto sabdw = sum_AB_dW_asym();
    auto sum_AB_dW = boost::yap::make_terminal(sabdw);

    store_old_value();
    auto ddts = ddt();
    auto ddto = boost::yap::make_terminal(ddts);

    solve(ddto(sum_AB_dW(-mp_ * (u_.b() - u_.a()))));

    // set iteration
    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, ParticleMass);
REGISTER_DEF_TYPE(TRANSPORTEQN, Conti);
REGISTER_DEF_TYPE(TRANSPORTEQN, TransientConti);
