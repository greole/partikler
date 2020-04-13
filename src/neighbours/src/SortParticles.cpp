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

#include "ObjectRegistry.hpp"

CountingSortParticles::CountingSortParticles(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg), pos_(objReg.get_pos()),
      sc_(objReg.get_object<Field<std::vector<SearchCube>>>("search_cubes")),
      si_(objReg.create_field<SizeTField>("sorting_idxs")),
      scd_(objReg.get_object<Generic<SearchCubeDomain>>("search_cube_domain")) {
}

void CountingSortParticles::execute() {

    log().info_begin() << "Sorting particles ";

    // PointField old_pos(pos_);

    // pos.store_old();

    // TODO too much copying
    auto [sc, si, pos] = countingSortParticles(scd_(), pos_);

    sc_ = sc;
    si_ = si;
    pos_ = pos;

    log().info_end();
    reorder_fields();
}

void CountingSortParticles::reorder_fields() {
    log().info_begin() << "Reordering particle fields ";

    for (auto &f : get_objReg().get_objects()) {
        if (f.second->get_name() == "Pos") continue;
        if (f.second->get_name() == "KernelW") continue;
        if (f.second->get_name() == "KerneldWdx") continue;
        if (f.first == "KerneldWdxNeighbour") continue;
        if (f.second->get_name() == "neighbour_pairs") continue;
        if (f.second->get_name() == "surface_dist") continue;
        if (f.second->get_name() == "search_cubes") continue;
        if (f.second->get_name() == "sorting_idxs") continue;
        if (f.second->get_name() == "idx") continue;

        auto type = f.second->get_type();
        std::cout << "reordering " << f.second->get_name() << std::endl;
        std::shared_ptr<SPHObject> *obj_ptr = &f.second;
        DISPATCH(obj_ptr, reorder_vector, type, si_);
    }

    log().info_end();
}

REGISTER_DEF_TYPE(SORTING, CountingSortParticles);
