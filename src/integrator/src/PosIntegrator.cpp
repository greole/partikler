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

#include "PosIntegrator.hpp"

#include "Time.hpp"

PosIntegrator::PosIntegrator(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : VectorFieldEquation(
          "Position", parameter, objReg, objReg.get_object<VectorField>("Pos")),
      id_(objReg.get_object<IntField>("id")),
      u_(objReg.get_object<VectorField>("u")),
      time_(objReg.get_object<TimeGraph>("TimeGraph")) {}

void PosIntegrator::execute() {

    log().info_begin() << " Updating Particle Positions";
    log().info() << " DeltaT " << time_.get_deltaT();

    // auto ddts = ddt();
    store_old_value();
    auto ddts = Ddt<VectorField>(time_.get_deltaT(), fo_, this->id_);
    auto ddto = boost::yap::make_terminal(ddts);

    solve(ddto(u_), false);

    // f_ += dx;

    // const float CFL = 0.01;
    // float current_max_dx = (pos_ - old_pos).norm().get_max();
    // // float dx_ratio = CFL*dx_max/current_max_dx;
    // float dx_ratio = CFL /current_max_dx;
    // float two = 2.0;
    // float change = min(two, dx_ratio);
    // time_.set_deltaT(min(time_.get_maxDeltaT(), time_.get_deltaT() *
    // change));

    log().info_end();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, PosIntegrator);
