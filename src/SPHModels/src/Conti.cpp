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

Conti::Conti (
    const std::string &model_name, YAML::Node parameter, RunTime &runTime):
    SPHModel(model_name, parameter, runTime),
    pos_(get_runTime().get_particle_positions()),
    np_(get_runTime().get_obj<SPHField<NeighbourPair>>("neighbour_pairs")),
    W_(get_runTime().get_obj<SPHFloatField>("KernelW")),
    rho_(get_runTime().get_obj<SPHFloatField>("rho")),
    lower_limit_(read_or_default_coeff<float>("lower_limit", 0.0))
{};

void Conti::execute() {

    log().info_begin() << "Computing density";

    rho_.set_uniform(0.0);

    rho_.weighted_sum(np_, W_);

    rho_.lower_limit(lower_limit_);

    log().info_end();
};

REGISTER_DEF_TYPE(TRANSPORTEQN, Conti);
