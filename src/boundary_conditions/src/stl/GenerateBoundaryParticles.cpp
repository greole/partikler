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

    // TODO  set boundaryIds:
    // 1. look in boundaryId map if boundary with given name exists
    //    if true get existing id
    //    if not calculate new id
    // 2. set boundary id to boundaryIds_ field
*/

#include "GenerateBoundaryParticles.hpp"

#include "Time.hpp"

GenerateBoundaryParticles::GenerateBoundaryParticles(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ParticleGeneratorBase(model_name, parameter, objReg),
      timeGraph_(local_objReg_.register_object<TimeGraph>(
          std::make_unique<TimeGraph>("TimeGraph", parameter, local_objReg_))),
      iterations_(read_coeff<int>("iterations")),
      write_freq_(read_or_default_coeff<int>("writeout", -1)),
      filename_(read_coeff<std::string>("file")),
      scale_(read_or_default_coeff<float>("scale", 1.0)),
      nu_(read_or_default_coeff<float>("nu", 1e-05)),
      rho_0_(read_or_default_coeff<float>("rho_0", 1.0)),
      p_0_(read_or_default_coeff<float>("p_0", 1e-05)),
      dxhratio_(read_or_default_coeff<float>("dxhratio", 1.05))
{}

YAML::Node GenerateBoundaryParticles::default_graph() {
    YAML::Node node; // starts out as NULL

    YAML::Node reader;
    // Reads the STL file
    reader["READER"]["model"] = "SPHSTLReader";
    reader["READER"]["file"] = filename_;

    node["pre"].push_back(reader);

    YAML::Node generator;
    // Generates the particles at the surface
    generator["GENERATOR"]["model"] = "SPHParticleGenerator";
    generator["GENERATOR"]["dx"] = dx_/scale_;
    generator["GENERATOR"]["id"] = fieldId_;
    node["pre"].push_back(generator);

    YAML::Node neighbours;
    // finds the neighbours
    neighbours["PARTICLENEIGHBOURS"]["model"] = "SPHSTLParticleNeighbours";
    neighbours["PARTICLENEIGHBOURS"]["dx"] = dx_/scale_ * dxhratio_;
    node["main"].push_back(neighbours);

    YAML::Node kernel;
    kernel["KERNEL"]["model"] = "STLWendland2D";
    kernel["KERNEL"]["h"] = dx_/scale_ * dxhratio_;
    node["main"].push_back(kernel);

    YAML::Node conti;
    conti["TRANSPORTEQN"]["model"] = "Conti";
    conti["TRANSPORTEQN"]["lower_limit"] = 0.001;
    node["main"].push_back(conti);

    YAML::Node pressure;
    pressure["TRANSPORTEQN"]["model"] = "Pressure";
    pressure["TRANSPORTEQN"]["rho_0"] = rho_0_ ;
    pressure["TRANSPORTEQN"]["p_0"] = p_0_;
    node["main"].push_back(pressure);

    YAML::Node visc;
    visc["TRANSPORTEQN"]["model"] = "Viscosity";
    visc["TRANSPORTEQN"]["nu"] = nu_;
    node["main"].push_back(visc);

    YAML::Node mom;
    mom["TRANSPORTEQN"]["model"] = "Momentum";
    node["main"].push_back(mom);

    YAML::Node integ;
    integ["TRANSPORTEQN"]["model"] = "STLPosIntegrator";
    node["main"].push_back(integ);

    if (debug()){
        YAML::Node writer;
        writer["EXPORT"]["model"] = "SuperSPHWriter";
        writer["EXPORT"]["name"] = name_;
        writer["EXPORT"]["writeout"] = write_freq_;
        node["main"].push_back(writer);
    }

    return node;
}


void GenerateBoundaryParticles::execute() {

    Logger logger {1};

    // Register Models

    auto node = default_graph();

    for (auto el : node["pre"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();

        auto model = ModelFactory::createInstance(
            model_namespace, model_name, model_name, it->second, local_objReg_);

        timeGraph_.push_back_pre(model);
    }

    timeGraph_.execute_pre();

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    // float kernel_relaxation = 1.0;
    // float noise_relaxation = 1.0;

    // surface slide particle have type 2
    // IntField &type = get_objReg().create_field<IntField>("type", 2);
    // SizeTField &idx = obj_reg.create_idx_field();

    // Set particle mass to match rho_0
    auto& spm = local_objReg_.get_object<Generic<float>>("specific_particle_mass");
    spm() = spm()*rho_0_;

    // Register Models
    for (auto el : node["main"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();

        auto model = ModelFactory::createInstance(
            model_namespace, model_name, model_name, it->second, local_objReg_);

        timeGraph_.push_back_main(model);
    }

    timeGraph_.execute_main();

    auto &loc_objs = local_objReg_.get_objects();

    auto& pos(local_objReg_.get_object<VectorField>("Pos"));

    scalePoints(pos, scale_);

    post_execute();
}

REGISTER_DEF_TYPE(BOUNDARY, GenerateBoundaryParticles);
