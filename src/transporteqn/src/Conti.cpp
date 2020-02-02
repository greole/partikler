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
          model_name,
          parameter,
          objReg,
          objReg.create_field<ScalarField>("rho", 0.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      lower_limit_(read_or_default_coeff<Scalar>("lower_limit", 0.0)),
      dx_(read_or_default_coeff<Scalar>("dx", 1.0)),
      mp_(
          objReg.get_object<Generic<Scalar>>("specific_particle_mass")()
          ) {

    log().info_begin() << "Updating specific particle mass";

    Scalar pi = (Scalar) M_PI;
    Scalar r= dx_/2.0;
    mp_ = rho_0_*4.0/3.0*r*r*r*pi;

    log().info_end();
}

void Conti::execute() {


    log().info_begin() << "Computing density";

    // TODO needs lazy reset of rho_;
    for (auto &el : f_) {
        el = 0;
    }

    m_sum_AB(mp_);
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

REGISTER_DEF_TYPE(TRANSPORTEQN, Conti);
