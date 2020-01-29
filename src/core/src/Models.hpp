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

#ifndef SPHMODEL_H
#define SPHMODEL_H

#include <algorithm>                // for max
#include <boost/hana/at.hpp>        // for at_t::operator()
#include <boost/yap/expression.hpp> // for as_expr, make_terminal
#include <iostream>                 // for operator<<, basic_ostream, endl
#include <map>                      // for map<>::iterator, operator!=, ope...
#include <memory>                   // for shared_ptr, allocator, __shared_...
#include <stdio.h>
#include <string>  // for string, operator<<, operator+
#include <utility> // for move, pair, make_pair
#include <vector>  // for vector

#include "Field.hpp"            // for ScalarField, KernelGradientField
#include "FieldOps.hpp"         // for sum_AB_impl
#include "Logger.hpp"           // for MSG, Logger
#include "Object.hpp"           // for SPHObject, ModelType
#include "ObjectRegistry.hpp"   // for ObjectRegistry
#include "SearchCubes.hpp"      // for NeighbourFieldAB
#include "yaml-cpp/node/impl.h" // for Node::~Node, Node::Node, Node::as
#include "yaml-cpp/node/node.h" // for Node
#include "yaml-cpp/yaml.h"

#define REGISTER_DEC_TYPE(NAME) static ModelRegister<NAME> reg

#define REGISTER_DEF_TYPE(CLASS, NAME)                                         \
    ModelRegister<NAME> NAME::reg(#CLASS "::" #NAME)

// Abstract base class for Models
class Model : public SPHObject {

  private:
    YAML::Node parameter_; // Parameter of given model

    ObjectRegistry &objReg_; // Reference to main runTime

    Logger log_; // Own logger instance for better scoping

    std::vector<std::shared_ptr<Model>> submodels_;

    bool debug_;

  public:
    Model(const std::string name, YAML::Node parameter, ObjectRegistry &objReg);

    YAML::Node get_parameter() { return parameter_; }

    template <class T> T read_coeff(std::string coeff_name) {
        T coeff = parameter_[coeff_name].as<T>();
        log().info() << "Reading coefficient " << coeff_name << " " << coeff;
        return coeff;
    }

    template <class T>
    T read_or_default_coeff(std::string coeff_name, T default_val) {
        if (parameter_[coeff_name]) return read_coeff<T>(coeff_name);
        log().info() << "Using default coefficient " << coeff_name << " "
                     << default_val;
        return default_val;
    }

    ObjectRegistry &get_objReg() { return objReg_; };

    Logger &log() { return logger_; };

    // TODO implement a pre and post execute member
    virtual void execute() = 0;

    void sub_model_push_back(std::shared_ptr<Model> m);

    void execute_submodels();

    bool debug() { return debug_; };
};

class ModelFactory {
    /**
       Factory for Models, implements a map of constructor pointers
       for all existing and registered models.
     */

    typedef std::map<
        std::string,
        Model *(*)(const std::string &, YAML::Node, ObjectRegistry &)>
        map_type;

  private:
    static map_type *map_;

  public:
    static std::shared_ptr<Model> createInstance(
        const std::string &model_type,
        const std::string &model_name,
        const std::string, //&objReg_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    static void print_models(const std::string &model_type);

  protected:
    static map_type *getMap() {
        // never delete'ed. (exist until program termination)
        // because we can't guarantee correct destruction order
        if (!map_) {
            map_ = new map_type;
        }
        return map_;
    }
};

template <typename T>
Model *createModel(
    // creates instance of Model at runTime
    const std::string &model_name, // model name at runtime e.g. wendland
    YAML::Node parameter,
    ObjectRegistry &objReg) {
    return new T(model_name, parameter, objReg);
}

template <typename T> struct ModelRegister : ModelFactory {

    ModelRegister(std::string const &s) {
        getMap()->insert(std::make_pair(s, &createModel<T>));
    }
};

// Default Models

class ModelGraph : public Model {

    // Defines a temporal order for submodels
    // submodels could also be realised via template args,
    // static ModelRegister<ModelGraph> reg;
    REGISTER_DEC_TYPE(ModelGraph);

  public:
    // TODO Move implementation to cpp file
    ModelGraph(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg)
        : Model(model_name, parameter, objReg) {};

    void execute() { execute_submodels(); };
};

// Base Class for particle generating models
//
// Base Class for particle generating models
// handles creation of local object registry
// and transferring
class ParticleGeneratorBase : public Model {

  protected:
    FieldIdMap &fieldIdMap_;

    ObjectRegistry local_objReg_;

    PointField &points_;

    VectorField &pos_;

    IntField &id_;

    std::string name_;

    int fieldId_;

    float dx_;

    Vec3 translation_vector_;

  public:
    ParticleGeneratorBase(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    Vec3 read_vector(YAML::Node parameter, std::string coeff);

    template <class T> void append(T &, std::string name) {
        auto &oreg = get_objReg();
        field_append(
            oreg.get_object<T>(name), local_objReg_.get_object<T>(name));
    }

    virtual void execute() {};

    void post_execute() {

        logger_.info_begin() << "Translate: " << translation_vector_;

        translatePoints(pos_, translation_vector_);

        logger_.info_end();

        // Transfer particles to main object registry
        // Here only the position and the boundary id and type id
        // need to be transferred
        logger_.info_begin() << "Transfering ";

        // TODO transfer the Vec3 version of Pos
        auto &oreg = get_objReg();
        if (oreg.object_exists("Pos")) {
            append(pos_, "Pos");
        } else {
            // Move the object if it doesn't exist in the main registry yet
            oreg.get_objects().push_back(*local_objReg_.get_object_ptr("Pos"));
        }

        id_.reserve(pos_.size());
        for (size_t i = 0; i < pos_.size(); i++) {
            id_.push_back(fieldId_);
        }

        if (oreg.object_exists("id")) {
            append(id_, "id");
        } else {
            // Move the object if it doesn't exist in the main registry yet
            oreg.get_objects().push_back(*local_objReg_.get_object_ptr("id"));
        }

        get_objReg().update_n_particles();
        logger_.info_end();
    };
};
#endif
