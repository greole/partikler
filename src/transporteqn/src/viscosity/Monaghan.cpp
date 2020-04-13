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
#include "Monaghan.hpp"

#include "Conti.hpp"
#include "Scalar.hpp"
#include "Time.hpp"

Monaghan::Monaghan(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorGradientEquation(
          "Viscosity",
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "tau",
              zero<VectorField::value_type>::val,
              {"taux", "tauy", "tauz"})),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      alpha_(read_or_default_coeff<Scalar>("alpha", 1e-05)),
      c_(read_or_default_coeff<Scalar>("c", 300)), u_(objReg.velocity()),
      mp_(objReg.get_object<ScalarField>("mp")), pos_(objReg.get_pos()) {

    // maxDt_ = 0.25 * h_ * h_ / nu_;
    // time_.set_model_timestep(model_name, maxDt_);
}

void Monaghan::execute() {

    log().info_begin() << "Computing dnu";

    // auto& u = momentum_.get();
    ScalarField &rho = conti_.get();

    // clang-format off
    auto &dW = this->get_objReg().template get_object<KernelGradientField>(
        "KerneldWdx");
    auto sum_AB_e = Sum_AB_asym<VectorField>(f_, np_);
    auto sum_AB = boost::yap::make_terminal(sum_AB_e);
    Scalar fact =  2.0*alpha_*h_*c_;
    // TODO do this via a solve method
    VectorFieldAB dist(this->np_.size(), {0,0,0});
    solve_inner_impl(this->np_, dist, pos_.a() - pos_.b());
    solve(
    sum_AB(
        // fact / (rho.a() + rho.b()) * ((ab(u_) * dist)/(dist*dist))),true
        mp_.b()*fact / (rho.a() + rho.b()) * ((dist*dW)/(dist*dist))*ab(u_)),true
        );
    // clang-format on

    log().info_end();
    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(VISCOSITY, Monaghan);
