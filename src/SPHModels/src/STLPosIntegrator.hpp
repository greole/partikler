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

#ifndef STLPOSINTEGRATOR_H
#define STLPOSINTEGRATOR_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"
#include "STLLimitedDx.hpp"

class STLPosIntegrator : public SPHModel {

    REGISTER_DEC_TYPE(STLPosIntegrator);

private:

    // const std::vector<Point> opoints_;
    const SPHIntField &type_;
    SPHField<Facet_handle> & facets_;
    const SPHSizeTField &idx_;

    // Out
    SPHVectorField &u_;
    SPHPointField& pos_;
    SPHGeneric<TimeInfo>& time_;

public:
    STLPosIntegrator(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
