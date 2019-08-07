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


#include "ParticleGenerator.hpp"

SPHSTLReader::SPHSTLReader(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)
    : SPHModel(model_name, parameter, runTime),
      fn_(parameter["file"].as<std::string>()) {


    log().info_begin() << "Reading input file: " <<  fn_;

    std::ifstream * istream= new std::ifstream(fn_);

    // read_STL(istream, points, facets, false);
    Polyhedron_builder_from_STL<HalfedgeDS> * builder = new
        Polyhedron_builder_from_STL<HalfedgeDS>  (*istream);
    log().info_end();

    log().info_begin() << "Constructing polyhedron";

    // TODO make
    // Create input polyhedron
    CGALPolyhedron *polyhedron = new CGALPolyhedron;
    polyhedron->delegate(*builder);

    get_runTime().get_obj_reg().register_object<SPHGeneric<CGALPolyhedron>>(
        std::make_unique<SPHGeneric<CGALPolyhedron>>(
            "polyhedron", "generic", *polyhedron
            )
        );

    log().info_end();
}

void SPHSTLReader::execute() {
    log().info() << "Doing nothing";
}

SPHParticleGenerator::SPHParticleGenerator(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)
    : SPHModel(model_name, parameter, runTime),
      polyhedron_(
          get_runTime().get_obj<SPHGeneric<CGALPolyhedron>>("polyhedron")),
      facets_(get_runTime().create_field<SPHField<Facet_handle>>("facets")),
      pos_(get_runTime().create_field<SPHPointField>("Pos")),

      dx_(parameter["dx"].as<float>()) {}

void SPHParticleGenerator::execute() {

    log().info_begin() << "Generating initial particles ";
    Generate_Points_at_Facets gpf(dx_, pos_.get_field(), facets_.get_field());

    std::for_each(
        polyhedron_().facets_begin(), polyhedron_().facets_end(), gpf);

    log().info_end() << "Generated " << pos_.size() << " Particles";

    get_runTime().update();
}

REGISTER_DEF_TYPE(READER, SPHSTLReader);
REGISTER_DEF_TYPE(GENERATOR, SPHParticleGenerator);
