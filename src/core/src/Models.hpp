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

#include "Datastructures.hpp"
#include "Object.hpp"
#include "ObjectRegistry.hpp"

#include "yaml-cpp/yaml.h"
#include <map>
#include <memory>
#include <stdio.h>

#define REGISTER_DEC_TYPE(NAME) static ModelRegister<NAME> reg

#define REGISTER_DEF_TYPE(CLASS, NAME) ModelRegister<NAME> NAME::reg(#CLASS "::" #NAME )

class Model : public SPHObject {
    /**
       Abstract base class for Models

     */

  private:
    YAML::Node parameter_; // Parameter of given model
    ObjectRegistry & objReg_;    // Reference to main runTime
    Logger log_;           // Own logger instance for better scoping

    std::vector<std::shared_ptr<Model>> submodels_;

  public:
    Model(
        const std::string name,
        YAML::Node parameter,
        ObjectRegistry & objReg)
        :
        SPHObject(name, "Model"),
        parameter_(parameter),
        objReg_(objReg),
        // TODO inherit verbosity from main Logger
        log_(Logger(3)),
        submodels_(
            std::vector<std::shared_ptr<Model>> ()
            )
    {
        logger_.set_scope(this->get_name());
        log().info() << "Creating "
              << this->get_type()
              << " "
              << this->get_name();
    }


    YAML::Node get_parameter() { return parameter_;}

    template<class T>
    T read_coeff(std::string coeff_name) {
        T coeff = parameter_[coeff_name].as<T>();
        log().info() << "Reading coefficient "
                     << coeff_name << " "
                     << coeff;
        return coeff;
    }

    template<class T>
    T read_or_default_coeff(std::string coeff_name, T default_val) {
        if (parameter_[coeff_name]) return read_coeff<T>(coeff_name);
        log().info() << "Using default coefficient "
                     << coeff_name << " "
                     << default_val;
        return default_val;
    }

    ObjectRegistry & get_objReg() {return objReg_;};

    Logger & log() {return logger_;};

    // TODO implement a pre and post execute member
    virtual void execute() = 0;

    void sub_model_push_back(std::shared_ptr<Model> m) {

        log().info()
            << " Registering: " << this->get_type() << " " << this->get_name();

        submodels_.push_back(m);
    };

    void execute_submodels () {
        log().info()
            << " Executing: " << this->get_type() << " " << this->get_name();

        for (auto model : submodels_) {
            log().info()
                      << " Executing: " << model->get_type() << " "
                      << model->get_name();
            model->execute();
        }
    }
};

struct ModelFactory {
    /**
       Factory for Models, implements a map of constructor pointers
       for all existing and registered models.
     */

    typedef std::map<
        std::string,
        Model *(*)(
            const std::string &,
            YAML::Node,
            ObjectRegistry&)> map_type;

  private:
    static map_type *map_;

  public:
    static std::shared_ptr<Model>
    createInstance(
        const std::string &model_type,
        const std::string &model_name,
        const std::string &objReg_name,
        YAML::Node parameter,
        ObjectRegistry &objReg) {
        const std::string delim = "::";
        map_type::iterator it = getMap()->find(model_type + delim + model_name);
        if (it == getMap()->end()) {
            std::cout << " no model named "
                      << model_name
                      << " found, available models for namespace "
                      << model_type
                      << std::endl;
            print_models(model_type);
            return 0;
        };
        return std::shared_ptr<Model>(
            it->second(objReg_name, parameter, objReg)
            );
    }

    static void print_models(
        const std::string &model_type
        ) {

        map_type::iterator it = getMap()->begin();
        while (it != getMap()->end())
        {
            std::string word = it->first;
            if (word.rfind(model_type, 0) == 0) {
                std::cout << word << std::endl;
            }

            it++;
        }
    }

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
    const std::string &model_name,  // model name at runtime e.g. wendland
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
        const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
        : Model(model_name, parameter, objReg) {};

    void execute() { execute_submodels(); };
};

// class TransportEqn : public Model {

// private:

//     NeighbourFieldAB& nb_;

// public: 

//     TransportEqn() {};

//     sum_AB() {};


// }

class TimeGraph : public Model {

    // Defines a temporal order for submodels
    // submodels could also be realised via template args,
    // static ModelRegister<ModelGraph> reg;
    REGISTER_DEC_TYPE(TimeGraph);

private:

    ModelGraph init_;

    ModelGraph main_;

    ModelGraph post_;

    int   current_timestep_;

    float current_time_;

    std::string name_;

    float endTime_;

    float deltaT_;

    float max_deltaT_;

    int iterations_;

    int write_freq_;

    int last_write_;

    bool iter_mode = false;


public:
    // TODO Move implementation to cpp file
    TimeGraph(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry & objReg)
        : Model(model_name, parameter, objReg),
          init_(ModelGraph("pre",  parameter, objReg)),
          main_(ModelGraph("main", parameter, objReg)),
          post_(ModelGraph("post", parameter, objReg)),
          current_timestep_(0),
          current_time_(0),
          name_(read_coeff<std::string>("name")),
          endTime_(read_or_default_coeff<float>("endTime", -1.0 )),
          deltaT_(read_or_default_coeff<float>("deltaT", -1.0 )),
          max_deltaT_(read_or_default_coeff<float>("max_deltaT", -1.0 )),
          iterations_(read_or_default_coeff<int>("iterations", 0)),
          write_freq_(read_or_default_coeff<int>("write_frequency", 0))
    {
        if (deltaT_ < 0) {iter_mode = true;}
    };

    void execute() {
        execute_pre();
        execute_main();
        execute_post();
    };

    void execute_pre() {
        init_.execute();
    }

    void execute_main() {
        // TODO register some kind of call back
        while (current_timestep_ < iterations_) {
            main_.execute();
            get_objReg().write_to_disk(current_timestep_, name_);
            current_timestep_++;
        };
    }

    void execute_post() {
        post_.execute();
    }


    void push_back_pre(std::shared_ptr<Model> m) {
        init_.sub_model_push_back(m);
    }

    void push_back_main(std::shared_ptr<Model> m) {
        main_.sub_model_push_back(m);
    }

    void push_back_post(std::shared_ptr<Model> m) {
        post_.sub_model_push_back(m);
    }

    float get_deltaT() {
        if (iter_mode) return 1.0e-32;
        return deltaT_;
    }

    float get_maxDeltaT() {
        if (iter_mode) return 1.0;
        return max_deltaT_;
    }

    void set_deltaT(float deltaT) {
        deltaT_ = deltaT;
    }

    int & get_current_timestep() {return current_timestep_;}
};


// class Ddt : Model {
//
//
// }



// TODO TransportedQuantity
// a class providing field to Equation depedency
// has an update method that calls the dependend transport eqn

// TODO TransportEqn Base class
// - Access previous timestep result if activated
// - has an update method
// - allows other TransportEqn as dependencies
// - Computes a d(...)

#endif
