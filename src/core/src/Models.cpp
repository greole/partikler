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

#include "Models.hpp"

Model::Model(
    const std::string name, YAML::Node parameter, ObjectRegistry &objReg)
    : SPHObject(name, ModelType), parameter_(parameter), objReg_(objReg),
      // TODO inherit verbosity from main Logger
      log_(Logger(3)), submodels_(std::vector<std::shared_ptr<Model>>()),
      debug_(read_or_default_coeff<bool>("debug", false)) {
    logger_.set_scope(this->get_name());
    log().info() << "Creating Model: " << this->get_name();
}

void Model::sub_model_push_back(std::shared_ptr<Model> m) {

    log().info() << "Registering: Submodel: " << this->get_name();

    submodels_.push_back(m);
}

void Model::execute_submodels() {
    log().info() << "Executing Submodel: " << this->get_type_str() << " "
                 << this->get_name();

    for (auto model : submodels_) {
        log().info() << "Executing: " << model->get_type_str() << " "
                     << model->get_name();
        model->execute();
    }
}

ModelFactory::map_type *ModelFactory::map_;

std::shared_ptr<Model> ModelFactory::createInstance(
    const std::string &model_type,
    const std::string &model_name,
    const std::string, //&objReg_name,
    YAML::Node parameter,
    ObjectRegistry &objReg) {
    std::cout << "createInstance" << model_name << std::endl;
    const std::string delim = "::";
    map_type::iterator it = getMap()->find(model_type + delim + model_name);
    if (it == getMap()->end()) {
        std::cout << " no model named " << model_name
                  << " found, available models for namespace " << model_type
                  << std::endl;
        print_models(model_type);
        // TODO throw if not found
    };

    Model *model_ptr = it->second(model_name, parameter, objReg);

    return objReg.register_object_get_ptr<Model>(
        std::shared_ptr<Model>(model_ptr));
}

void ModelFactory::print_models(const std::string &model_type) {

    map_type::iterator it = getMap()->begin();
    while (it != getMap()->end()) {
        std::string word = it->first;
        if (word.rfind(model_type, 0) == 0) {
            std::cout << word << std::endl;
        }

        it++;
    }
}

ParticleGeneratorBase::ParticleGeneratorBase(
    const std::string &model_name,
    YAML::Node parameter,
    ObjectRegistry &objReg):
    Model(model_name, parameter, objReg),
    fieldIdMap_(objReg.get_object<FieldIdMap>("FieldIdMap")),
    local_objReg_(ObjectRegistry()),
    points_(local_objReg_.create_field<PointField>("Points", {}, {"X", "Y", "Z"})),
    pos_(local_objReg_.create_field<VectorField>("Pos", {}, {"X", "Y", "Z"})),
    id_(local_objReg_.create_field<IntField>("id")),
    name_(read_coeff<std::string>("name")),
    fieldId_(fieldIdMap_.append(name_)),
    dx_(read_coeff<float>("dx")),
    translation_vector_(read_vector(parameter, "translate"))
{}

Vec3
ParticleGeneratorBase::read_vector(YAML::Node parameter, std::string coeff) {

    if (!parameter[coeff]) return {0, 0, 0};

    auto p = parameter[coeff];
    return {p[0].as<float>(), p[1].as<float>(), p[2].as<float>()};
}

REGISTER_DEF_TYPE(CORE, ModelGraph);
