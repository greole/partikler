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

#include "AkinciCurvature.hpp"

#include "Conti.hpp"
#include "Time.hpp"
#include <math.h>

AkinciCurvature::AkinciCurvature(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          "Curvature",
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "fCurv", {}, {"fCurvx", "fCurvy", "fCurvz"})),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      mp_(objReg.get_object<ScalarField>("mp")),
      gamma_(read_or_default_coeff<Scalar>("gamma", 1.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      n_(objReg.get_object<VectorField>("normal")) {}

void AkinciCurvature::execute() {

    log().info_begin() << "Computing coehesive forces";
    iteration_ = time_.get_current_timestep();

    Scalar gamma = -1.0 * gamma_;

    ScalarField &rho = conti_.get();
    auto saba = sum_AB_a();
    auto sum_AB_e = boost::yap::make_terminal(saba);

    solve(
        sum_AB_e(2.0 * rho_0_ / (rho.a() + rho.b()) * (mp_.a() * gamma * ab(n_))));

    log().info_end();

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(SURFACETENSION, AkinciCurvature);
