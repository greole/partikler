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

#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"


class Viscosity : public SPHModel {

    REGISTER_DEC_TYPE(Viscosity);

private:

    // In
    const SPHField<searchcubes::NeighbourPair> &np_;
    const SPHField<VectorPair> &dW_;
    const SPHVectorField &u_;
    const SPHPointField &pos_; // Particle positions

    // Out
    SPHVectorField &dnu_;

    // Coeffs
    float nu_;

public:
    Viscosity(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif