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
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      dx_(parameter["dx"].as<float>()),
      pos_(objReg.get_particle_positions()),
      facets_(objReg.get_object<Field<Facet_handle>>("facets")),
      sc_(objReg.create_field<Field<searchcubes::SearchCube>>(
          "search_cubes")),
      np_(objReg.create_field<Field<searchcubes::NeighbourPair>>(
          "neighbour_pairs")),
      sd_(objReg.create_field<Field<STLSurfaceDist>>("surface_dist")),
      scd_(objReg
               .create_generic<Generic<searchcubes::SearchCubeDomain>>(
                   "search_cube_domain")) {

    auto pcs = ModelFactory::createInstance(
        "SORTING", "CountingSortParticles", "sort", parameter, objReg);

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
