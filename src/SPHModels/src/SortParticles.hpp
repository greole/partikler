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

#ifndef SORTPARTRICLES_H
#define SORTPARTRICLES_H

#include "SPHDatastructures.hpp"
#include "SPHModels.hpp"

#include "SearchCubes.hpp"

class CountingSortParticles : public SPHModel {

    REGISTER_DEC_TYPE(CountingSortParticles);

  private:
    // In
    SPHPointField &pos_;

    SPHField<searchcubes::SearchCube> &sc_;

    //  Sorting indexes
    SPHSizeTField &si_;

    SPHGeneric<searchcubes::SearchCubeDomain> &scd_;

  public:
    CountingSortParticles(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();

    void reorder_fields();
};

#endif