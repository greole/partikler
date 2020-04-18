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

#include "Field.hpp"    // for ScalarField, KernelGradientField
#include "FieldOps.hpp" // for sum_AB_impl
#include "Logger.hpp"   // for MSG, Logger
#include "Object.hpp"   // for SPHObject, ModelType
#include "yaml-cpp/yaml.h"

#include "Scalar.hpp"
#include "Vec3.hpp"
#include "static_block.hpp"

class ObjectRegistry;

#define REGISTER_DEC_TYPE(NAME)                                                \
  public:                                                                      \
    static ModelRegister<NAME> reg

#define REGISTER_DEF_TYPE(CLASS, NAME)                                         \
    ModelRegister<NAME> NAME::reg(#CLASS "::" #NAME)

template <class T>
T read_coeff_impl(YAML::Node const &parameter, std::string coeff_name) {
    T coeff = parameter[coeff_name].as<T>();
    return coeff;
}

template <class T>
T read_or_default_coeff_impl(
    YAML::Node const &parameter, std::string coeff_name, T default_val) {
    if (parameter[coeff_name]) return read_coeff_impl<T>(parameter, coeff_name);
    return default_val;
}

// Abstract base class for Models
class Model : public SPHObject {

  protected:
    YAML::Node parameter_; // Parameter of given model

    ObjectRegistry &objReg_; // Reference to main runTime

    Logger log_; // Own logger instance for better scoping

    std::vector<std::shared_ptr<Model>> submodels_;

    bool debug_;

  public:
    Model(const std::string name, YAML::Node parameter, ObjectRegistry &objReg);

    YAML::Node get_parameter() { return parameter_; }

    template <class T> T read_coeff(std::string coeff_name) {
        T coeff = read_coeff_impl<T>(parameter_, coeff_name);
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

    Vec3 read_vec3(std::string coeff_name) {
        Scalar x = parameter_[coeff_name][0].as<Scalar>();
        Scalar y = parameter_[coeff_name][1].as<Scalar>();
        Scalar z = parameter_[coeff_name][2].as<Scalar>();
        Vec3 ret {x, y, z};
        log().info() << "Reading coefficient " << coeff_name << " " << ret;
        return ret;
    }

    Vec3 read_or_default_vec3(std::string coeff_name, Vec3 default_val) {
        if (parameter_[coeff_name]) return read_vec3(coeff_name);
        log().info() << "Using default coefficient " << coeff_name << " "
                     << default_val;
        return default_val;
    }

    ObjectRegistry &get_objReg() { return objReg_; };

    Logger &log() { return logger_; };

    // TODO implement a pre and post execute member
    virtual void execute() = 0;

    void sub_model_push_back(std::shared_ptr<Model> m);

    std::vector<std::shared_ptr<Model>> &get_submodels() { return submodels_; }

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
#endif
