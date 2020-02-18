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

#include "Scalar.hpp"
#include "Time.hpp"
#include "Conti.hpp"


Viscosity::Viscosity(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorGradientEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "tau",
              zero<VectorField::value_type>::val,
              {"taux", "tauy", "tauz"})),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      nu_(read_or_default_coeff<Scalar>("nu", 1e-05)),
      u_(objReg.velocity()),
      mp_(objReg.get_object<Generic<Scalar>>("specific_particle_mass")()),
      pos_(objReg.get_pos()) {

    maxDt_ = 0.25 * h_*h_/nu_;
    time_.set_model_timestep(model_name, maxDt_);
}

void Viscosity::execute() {

    log().info_begin() << "Computing dnu";

    // auto& u = momentum_.get();
    int it = time_.get_current_timestep();
    ScalarField &rho = conti_.get(it);

    // clang-format off
    auto normSqr = boost::yap::make_terminal(NormSqr_Wrapper());
    Scalar fact =  mp_*10.0*nu_;
    // sum_AB_dW(
    //     fact * rho * (ab(u_) * ab(pos_))
    //   / ( ( rho.a() + rho.b() ) * ( normSqr(ab(pos_)) ) )
    //     );
    // clang-format on

    log().info_end();
    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Viscosity);
