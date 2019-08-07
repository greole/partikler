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

#include "ParticleNeighbours.hpp"

SPHSTLParticleNeighbours::SPHSTLParticleNeighbours(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)
    : SPHModel(model_name, parameter, runTime),
      dx_(parameter["dx"].as<float>()),
      pos_(get_runTime().get_particle_positions()),
      facets_(get_runTime().get_obj<SPHField<Facet_handle>>("facets")),
      sc_(get_runTime().create_field<SPHField<searchcubes::SearchCube>>(
          "search_cubes")),
      np_(get_runTime().create_field<SPHField<searchcubes::NeighbourPair>>(
          "neighbour_pairs")),
      sd_(get_runTime().create_field<SPHField<STLSurfaceDist>>("surface_dist")),
      scd_(get_runTime()
               .create_generic<SPHGeneric<searchcubes::SearchCubeDomain>>(
                   "search_cube_domain")) {

    auto pcs = SPHModelFactory::createInstance(
        "SORTING", "CountingSortParticles", "sort", parameter, runTime);

    sub_model_push_back(pcs);
};

void SPHSTLParticleNeighbours::execute() {

    // log().set_scope("")
    log().info() << " Constructing particle neighbours";
    update_search_cube_domain();

    execute_submodels();

    // TODO move sorting to countingSortParticles submodel
    // reorder_vector(
    //     get_runTime().get_obj<SPHSizeTField>("sorting_idxs").get_field(),
    //     facets_.get_field());

    log().info_begin() << "Call createSTLNeighbours";

    // TODO to much copying
    auto [np, sd] = createSTLNeighbours(
        scd_(),
        pos_.get_field(),
        sc_.get_field(),
        facets_.get_field());

    np_.set_field(np);
    sd_.set_field(sd);

    log().info() << "Found " << np_.size() << " Neighbour Pairs";

    log().info_end() << "DONE";
};

void SPHSTLParticleNeighbours::update_search_cube_domain() {
    scd_() = searchcubes::initSearchCubeDomain(
        pos_.get_field(), search_cube_size_*dx_);
}

REGISTER_DEF_TYPE(PARTICLENEIGHBOURS, SPHSTLParticleNeighbours);
