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

#include "SortParticles.hpp"

CountingSortParticles::CountingSortParticles(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)
    : SPHModel(model_name, parameter, runTime),
      pos_(get_runTime().get_particle_positions()),
      sc_(get_runTime().get_obj<SPHField<searchcubes::SearchCube>>(
          "search_cubes")),
      si_(get_runTime().create_field<SPHSizeTField>("sorting_idxs")),
      scd_(get_runTime().get_obj<SPHGeneric<searchcubes::SearchCubeDomain>>(
          "search_cube_domain")) {};

void CountingSortParticles::execute(){

    log().info_begin() << "Sorting particles ";

    std::vector<Point> old_pos = pos_.get_field();

    // TODO too much copying
    auto [sc, si, pos] = countingSortParticles(
        scd_(),
        old_pos
        );

    sc_.set_field(sc);
    si_.set_field(si);
    pos_.set_field(pos);

    log().info_end();
    reorder_fields();
};

void CountingSortParticles::reorder_fields(){
    log().info_begin() << "Reordering particle fields ";

    for (auto &f : get_runTime().get_obj_reg().get_objects()) {
        if (f->get_name() == "Pos" ) continue;
        if (f->get_name() == "KernelW" ) continue;
        if (f->get_name() == "KerneldWdx" ) continue;
        if (f->get_name() == "neighbour_pairs" ) continue;
        if (f->get_name() == "surface_dist" ) continue;
        if (f->get_name() == "search_cubes" ) continue;
        if (f->get_name() == "sorting_idxs" ) continue;
        f->reorder(si_.get_field());
    }

    log().info_end();
}

REGISTER_DEF_TYPE(SORTING, CountingSortParticles);
