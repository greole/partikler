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

# include "STLPosIntegrator.hpp"


STLPosIntegrator::STLPosIntegrator(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime),
      type_(get_runTime().get_obj<SPHIntField>("type")),
      facets_(get_runTime().get_obj<SPHField<Facet_handle>>("facets")),
      idx_(get_runTime().get_obj<SPHSizeTField>("idx")),
      u_(get_runTime().get_obj<SPHVectorField>("u")),
      pos_(get_runTime().get_obj<SPHPointField>("Pos")),
      time_(get_runTime().get_obj<SPHGeneric<TimeInfo>>("TimeInfo"))
{};

void STLPosIntegrator::execute() {

    log().info_begin() << " Updating Particle Positions";
    log().info() << " DeltaT " << time_().deltaT;

    SPHPointField old_pos(pos_);

    SPHVectorField dx = STL_limited_dx(u_, time_().deltaT, facets_, type_, idx_, pos_);

    pos_ += dx;

    const float CFL = 0.01;
    float current_max_dx = (pos_ - old_pos).norm().get_max();
    // float dx_ratio = CFL*dx_max/current_max_dx;
    float dx_ratio = CFL /current_max_dx;
    float two = 2.0;
    float change = min(two, dx_ratio);
    time_().deltaT =  min(time_().maxDeltaT, time_().deltaT * change);

    log().info_end();
};

REGISTER_DEF_TYPE(TRANSPORTEQN, STLPosIntegrator);
