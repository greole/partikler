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



#include "GenerateBoundaryParticles.hpp"

GenerateBoundaryParticles::GenerateBoundaryParticles(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      local_objReg_(ObjectRegistry()),
      timeGraph_(local_objReg_.register_object<TimeGraph>(
                     std::make_unique<TimeGraph>("TimeGraph", parameter, local_objReg_))),
      boundaryIds_(local_objReg_.create_field<IntField>("boundary")),
      typeIds_(local_objReg_.create_field<IntField>("type")),
      idx_(local_objReg_.create_field<SizeTField>("idx")),
      pos_(local_objReg_.create_field<PointField>("Pos")),
      iterations_(read_coeff<int>("iterations")),
      write_freq_(read_or_default_coeff<int>("writeout", -1)),
      filename_(read_coeff<std::string>("file")),
      boundary_name_(read_coeff<std::string>("name")),
      dx_(read_coeff<float>("dx"))
{};

YAML::Node GenerateBoundaryParticles::default_graph() {
    YAML::Node node;  // starts out as NULL

    YAML::Node reader;
    reader["READER"]["model"] = "SPHSTLReader";
    reader["READER"]["file"] = filename_;

    node["pre"].push_back(reader);

    YAML::Node generator;
    generator["GENERATOR"]["model"] = "SPHParticleGenerator";
    generator["GENERATOR"]["dx"] = dx_;
    node["pre"].push_back(generator);

    YAML::Node neighbours;
    neighbours["PARTICLENEIGHBOURS"]["model"] = "SPHSTLParticleNeighbours";
    neighbours["PARTICLENEIGHBOURS"]["dx"] = dx_ * 1.05;
    node["main"].push_back(neighbours);

    YAML::Node kernel;
    kernel["KERNEL"]["model"] = "STLWendland2D";
    kernel["KERNEL"]["h"] = dx_ * 1.05;
    node["main"].push_back(kernel);

    YAML::Node conti;
    conti["TRANSPORTEQN"]["model"] = "Conti";
    conti["TRANSPORTEQN"]["lower_limit"] = 0.001;
    node["main"].push_back(conti);

    YAML::Node pressure;
    pressure["TRANSPORTEQN"]["model"] = "Pressure";
    node["main"].push_back(pressure);

    YAML::Node visc;
    visc["TRANSPORTEQN"]["model"] = "Viscosity";
    visc["TRANSPORTEQN"]["nu"] = 100;
    node["main"].push_back(visc);

    YAML::Node mom;
    mom["TRANSPORTEQN"]["model"] = "Momentum";
    node["main"].push_back(mom);

    YAML::Node integ;
    integ["TRANSPORTEQN"]["model"] = "STLPosIntegrator";
    node["main"].push_back(integ);

    return node;
}

void GenerateBoundaryParticles::execute() {

    Logger logger {1};

    // Register Models

    auto node = default_graph();

    for (auto el: node["pre"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();

        auto model = ModelFactory::createInstance(
            model_namespace,
            model_name,
            model_name,
            it->second,
            local_objReg_
            );

        timeGraph_.push_back_pre(model);
    }

    timeGraph_.execute_pre();

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    // float kernel_relaxation = 1.0;
    float noise_relaxation = 1.0;


    // surface slide particle have type 2
    // IntField &type = get_objReg().create_field<IntField>("type", 2);
    // SizeTField &idx = obj_reg.create_idx_field();

    // Register Models
    for (auto el: node["main"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();

        auto model = ModelFactory::createInstance(
            model_namespace,
            model_name,
            model_name,
            it->second,
            local_objReg_
            );

        timeGraph_.push_back_main(model);
    }

    timeGraph_.execute_main();

    auto& loc_objs = local_objReg_.get_objects();

    logger_.info_begin() << "Transfering ";

    // check if  field already present
    for (auto &obj: loc_objs) {
        auto name = obj->get_name();
        auto type = obj->get_type();
        // "int" IntField
        // "long" SizeTField
        // "float" FloatField
        // "Point" PointField
        // "vector" VectorField
        // TODO Refactor this

        std::cout << "Transfering" << name << std::endl;

        if (get_objReg().object_exists(name)) {

            if (type=="int") {

                get_objReg().get_object<IntField>(name).append(
                    local_objReg_.get_object<IntField>(name).get_field());
            }

            if (type=="long") {

                get_objReg().get_object<SizeTField>(name).append(
                    local_objReg_.get_object<SizeTField>(name).get_field());
            }

            if (type=="float") {

                get_objReg().get_object<FloatField>(name).append(
                    local_objReg_.get_object<FloatField>(name).get_field());
            }

            if (type=="vector") {

                get_objReg().get_object<VectorField>(name).append(
                    local_objReg_.get_object<VectorField>(name).get_field());
            }

            if (type=="Point") {

                get_objReg().get_object<PointField>(name).append(
                    local_objReg_.get_object<PointField>(name).get_field());
            }

        } else {
            get_objReg().get_objects().push_back(
                std::move(obj)
                );
        }
    }

    logger_.info_end();

    // write
};

REGISTER_DEF_TYPE(BOUNDARY, GenerateBoundaryParticles);
