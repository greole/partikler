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

#include "BrownianMotion.hpp"

BrownianMotion::BrownianMotion(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime),
      u_(get_runTime().get_obj<SPHVectorField>("u")),
      dx_(read_or_default_coeff<float>("dx", 1.0)) {};

void BrownianMotion::execute() {
    srand(time(NULL));
    FOR_ALL(u_, ctr) {
        float x = ((float)(rand() % 100)) / 100.0 - 0.5;
        float y = ((float)(rand() % 100)) / 100.0 - 0.5;
        float z = ((float)(rand() % 100)) / 100.0 - 0.5;

        Vector uv = u_[ctr];

        Vector rv {uv[0] + x * dx_, uv[1] + y * dx_, uv[2] + z * dx_};
        u_[ctr] = rv;
    }
};

REGISTER_DEF_TYPE(TRANSPORTEQN, BrownianMotion);
