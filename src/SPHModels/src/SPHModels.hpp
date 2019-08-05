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

#include "SPHDatastructures.hpp"
#include "RunTime.hpp"

#include "yaml-cpp/yaml.h"
#include <map>
#include <stdio.h>

#define REGISTER_DEC_TYPE(NAME) static ModelRegister<NAME> reg

#define REGISTER_DEF_TYPE(CLASS, NAME) ModelRegister<NAME> NAME::reg(#CLASS "::" #NAME )

class SPHModel : public SPHObject {
    /**
       Abstract base class for SPHModels

     */

  private:
    YAML::Node parameter_; // Parameter of given model
    RunTime & runTime_;    // Reference to main runTime
    Logger log_;           // Own logger instance for better scoping

    std::vector<SPHModel *> submodels_ = {};

  public:
    SPHModel(
        const std::string name,
        YAML::Node parameter,
        RunTime & runTime)
        :
        SPHObject(name, "Model",  true),
        parameter_(parameter),
        runTime_(runTime),
        // TODO inherit verbosity from main Logger
        log_(Logger(3)) {
        std::cout << " creating "
                  << this->get_type()
                  << " "
                  << this->get_name()
                  << std::endl;
        logger_.set_scope(this->get_name());
        // std::cout << parameter_.as<std::string>() << std::endl;
    };

    template<class T>
    T read_coeff(std::string coeff_name) {
        T coeff = parameter_[coeff_name].as<T>();
        log().info() << "Reading coefficient "
                     << coeff_name << " "
                     << coeff;
        return coeff;
    };

    template<class T>
    T read_or_default_coeff(std::string coeff_name, T default_val) {
        if (parameter_[coeff_name]) return read_coeff<T>(coeff_name);
        log().info() << "Using default coefficient "
                     << coeff_name << " "
                     << default_val;
        return default_val;
    };

    RunTime & get_runTime() {return runTime_;};

    Logger & log() {return logger_;};

    // TODO implement a pre and post execute member
    virtual void execute() = 0;

    void sub_model_push_back(SPHModel &m) {

        std::cout << " Registering: " << this->get_type() << " "
                  << this->get_name() << std::endl;

        submodels_.push_back(&m);
    };


    void execute_submodels () {
        std::cout << " Executing: " << this->get_type() << " "
                  << this->get_name() << std::endl;

        for (auto model : submodels_) model->execute();
    }
};

struct SPHModelFactory {
    /**
       Factory for SPHModels, implements a map of constructor pointers
       for all existing and registered models.
     */

    typedef std::map<
        std::string,
        SPHModel *(*)(
            const std::string &,
            YAML::Node,
            RunTime&)> map_type;

  private:
    static map_type *map_;

  public:
    static SPHModel *
    createInstance(
        const std::string &model_type,
        const std::string &model_name,
        const std::string &runTime_name,
        YAML::Node parameter,
        RunTime &runTime) {
        const std::string delim = "::";
        map_type::iterator it = getMap()->find(model_type + delim + model_name);
        if (it == getMap()->end()) {
            std::cout << " no model named "
                      << model_name
                      << " found, available models:"
                      << std::endl;
            print_models(model_type);
            return 0;
        };
        return it->second(runTime_name, parameter, runTime);
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
SPHModel *createModel(
    // creates instance of SPHModel at runTime
    const std::string &model_name,  // model name at runtime e.g. wendland
    YAML::Node parameter,
    RunTime &runTime) {
    return new T(model_name, parameter, runTime);
}

template <typename T> struct ModelRegister : SPHModelFactory {

    ModelRegister(std::string const &s) {
        getMap()->insert(std::make_pair(s, &createModel<T>));
    }
};

// Default Models

class SPHModelGraph : public SPHModel {

    // Defines a temporal order for submodels
    // submodels could also be realised via template args,
    // static ModelRegister<SPHModelGraph> reg;
    REGISTER_DEC_TYPE(SPHModelGraph);

  public:
    // TODO Move implementation to cpp file
    SPHModelGraph(
        const std::string &model_name,
        YAML::Node parameter,
        RunTime & runTime)
        : SPHModel(model_name, parameter, runTime) {};

    void execute() {
        execute_submodels();
    };
};

#endif
