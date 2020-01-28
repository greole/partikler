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
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      fn_(parameter["file"].as<std::string>()) {

    log().info_begin() << "Reading input file: " << fn_;

    std::ifstream *istream = new std::ifstream(fn_);

    // read_STL(istream, points, facets, false);
    Polyhedron_builder_from_STL<HalfedgeDS> *builder =
        new Polyhedron_builder_from_STL<HalfedgeDS>(*istream);

    // objReg.register_object<Generic<Polyhedron_builder_from_STL<HalfedgeDS>>>(
    //     std::make_unique<Generic<Polyhedron_builder_from_STL<HalfedgeDS>>>(
    //         "builder", "generic", *builder
    //         )
    //     );

    log().info_end();

    log().info_begin() << "Constructing polyhedron";

    // TODO make
    // Create input polyhedron
    CGALPolyhedron *polyhedron = new CGALPolyhedron;

    polyhedron->delegate(*builder);

    objReg.register_object<Generic<CGALPolyhedron>>(
        std::make_unique<Generic<CGALPolyhedron>>(
            "polyhedron", GenericType, *polyhedron));

    delete (polyhedron);
    delete (builder);
    delete (istream);

    log().info_end();
}

void SPHSTLReader::execute() { log().info() << "Doing nothing"; }

SPHParticleGenerator::SPHParticleGenerator(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg), boundary_id_(read_coeff<int>("id")),
      polyhedron_(objReg.get_object<Generic<CGALPolyhedron>>("polyhedron")),
      facets_(objReg.create_field<Field<std::vector<Facet_handle>>>("facets")),
      points_(objReg.create_field<PointField>("Points", {}, {"X", "Y", "Z"})),
      pos_(objReg.create_field<VectorField>("Pos", {}, {"X", "Y", "Z"})),
      idx_(objReg.create_field<SizeTField>("idx")),
      type_(objReg.create_field<IntField>("type")),
      boundary_(objReg.create_field<IntField>("boundary")),
      dx_(parameter["dx"].as<float>()) {}

void SPHParticleGenerator::execute() {

    log().info_begin() << "Generating initial particles ";
    Generate_Points_at_Facets gpf(dx_, points_, facets_);

    size_t n_0 = points_.size();

    std::for_each(
        polyhedron_().facets_begin(), polyhedron_().facets_end(), gpf);

    size_t n = points_.size();

    // TODO find a better solution
    for (auto &p : points_) {
        pos_.push_back({(float)p[0], (float)p[1], (float)p[2]});
    }

    log().info_end() << "Generated " << n << " Particles";

    get_objReg().set_n_particles(n);

    // TODO START REFACTOR THIS
    Compute_Facet_Area ca;

    log().info_begin() << "Computing Surface Area";

    float surface_area = (float)std::accumulate(
        boost::make_transform_iterator(polyhedron_().facets_begin(), ca),
        boost::make_transform_iterator(polyhedron_().facets_end(), ca),
        0.);

    float volume = surface_area * dx_;
    log().info_end() << "Total Surface Area" << surface_area;

    get_objReg().register_object<Generic<float>>(
        std::make_unique<Generic<float>>(
            "specific_particle_mass", GenericType, volume / n));

    for (size_t i = 0; i < (n - n_0); i++) {
        idx_.push_back(n_0 + i);
        type_.push_back(2);
        boundary_.push_back(boundary_id_);
    }
}

REGISTER_DEF_TYPE(READER, SPHSTLReader);
REGISTER_DEF_TYPE(GENERATOR, SPHParticleGenerator);
