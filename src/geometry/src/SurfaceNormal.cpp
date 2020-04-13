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

#include "SurfaceNormal.hpp"

#include "Time.hpp"

SurfaceNormal::SurfaceNormal(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          "SurfaceNormal",
          parameter,
          objReg,
          objReg.create_field<VectorField>("normal", {}, {"nx", "ny", "nz"})),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      mp_(objReg.get_object<ScalarField>("mp")),
      h_(read_or_default_coeff<Scalar>("h", 1.0)) {}

void SurfaceNormal::execute() {

    ScalarField &rho = conti_.get();
    auto sum_AB_dW = boost::yap::make_terminal(sum_AB_dW_s);

    solve(h_ * sum_AB_dW(mp_.b() / rho.b()));

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(SURFACETENSION, SurfaceNormal);
