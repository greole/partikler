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

#include "STLParticleNeighbours.hpp"

#include "ObjectRegistry.hpp"
#include "cgal/CGALSearchCubes.hpp" // for SearchCube, NeighbourFieldAB, Searc...

SPHSTLParticleNeighbours::SPHSTLParticleNeighbours(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg), dx_(parameter["dx"].as<float>()),
      points_(objReg.get_object<PointField>("points")),
      facets_(objReg.get_object<Field<std::vector<Facet_handle>>>("facets")),
      sc_(objReg.create_field<SearchCubeFieldAB>("search_cubes")),
      np_(objReg.create_field<NeighbourFieldAB>("neighbour_pairs")),
      sd_(objReg.create_field<STLSurfaceDistAB>("surface_dist")),
      scd_(objReg.create_generic<Generic<SearchCubeDomain>>(
          "search_cube_domain")) {

    auto pcs = ModelFactory::createInstance(
        "SORTING", "CountingSortParticles", "sort", parameter, objReg);

    sub_model_push_back(pcs);
}

void SPHSTLParticleNeighbours::execute() {

    // log().set_scope("")
    log().info() << " Constructing particle neighbours";
    update_search_cube_domain();
    execute_submodels();

    // TODO move sorting to countingSortParticles submodel
    reorder_vector(
        facets_, get_objReg().get_object<SizeTField>("sorting_idxs"));

    log().info_begin() << "Call createSTLNeighbours";

    // TODO to much copying
    auto [np, sd] = createSTLNeighbours(scd_(), points_, sc_, facets_);

    // np_.set_field(np);
    np_ = np;
    // sd_.set_field(sd);
    sd_ = sd;

    log().info() << "Found " << np_.size() << " Neighbour Pairs";

    log().info_end() << "DONE";
}

void SPHSTLParticleNeighbours::update_search_cube_domain() {
    scd_() = initSearchCubeDomain(
        cgal_bounding_box(points_), search_cube_size_ * dx_);
}

REGISTER_DEF_TYPE(PARTICLENEIGHBOURS, SPHSTLParticleNeighbours);
