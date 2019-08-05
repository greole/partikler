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


#ifndef Pressure_H
#define Pressure_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"

class Pressure : public SPHModel {

    REGISTER_DEC_TYPE(Pressure);

  private:

    // In
    // Density
    const SPHFloatField &rho_;
    const SPHField<searchcubes::NeighbourPair> &np_;
    const SPHFloatField &W_;
    const SPHField<VectorPair> &dW_;

    // Out
    // Pressure
    SPHFloatField &p_;
    SPHVectorField &dp_;

    // Coeffs
    const float c_;
    const float rho_0_;
    const float gamma_;
    const float p_0_;
    const float prefac_;

  public:

    Pressure(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
