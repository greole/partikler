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

Conti::Conti(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ScalarFieldEquation(
          "Conti",
          parameter,
          objReg,
          objReg.create_field<ScalarField>("rho", 0.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      lower_limit_(read_or_default_coeff<Scalar>("lower_limit", 0.0)),
      mp_(objReg.get_object<Generic<Scalar>>("specific_particle_mass")()) {}

void Conti::execute() {

    log().info_begin() << "Computing density";

    auto sum_AB_e = boost::yap::make_terminal(sum_AB_s);

    solve(
        mp_*sum_AB_e(W_),
        true
        );

    // TODO do it lazily
    for (auto &el : f_) {
        if (el < lower_limit_) {
            el = lower_limit_;
        }
    }

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
      u_(objReg.velocity()),
      mp_(objReg.get_object<Generic<Scalar>>("specific_particle_mass")()) {}

void TransientConti::execute() {

    log().info_begin() << "Computing density";

    auto &dW = get_objReg().template get_object<KernelGradientField>("KerneldWdx");

    auto sum_AB_s = Sum_AB_dW_sym<ScalarField>(f_, np_, dW, dW);
    auto sum_AB_e = boost::yap::make_terminal(sum_AB_s);

    auto ddt = boost::yap::make_terminal(ddt_);

    solve(
        ddt(sum_AB_e(-mp_*(u_.b()-u_.a()))),
        true
        );

    // set iteration
    iteration_ = time_.get_current_timestep();

    log().info_end();
}


REGISTER_DEF_TYPE(TRANSPORTEQN, Conti);
REGISTER_DEF_TYPE(TRANSPORTEQN, TransientConti);
