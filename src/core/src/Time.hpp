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

#ifndef PARTIKLER_TIME_INCLUDED
#define PARTIKLER_TIME_INCLUDED

#include "Models.hpp"
#include "ObjectRegistry.hpp" // for ObjectRegistry
#include "Scalar.hpp"
#include "yaml-cpp/yaml.h"

// TODO use std::variant<float,size_t> to
class TimeGraph : public Model {

    // Defines a temporal order for submodels
    // submodels could also be realised via template args,
    // static ModelRegister<ModelGraph> reg;
    REGISTER_DEC_TYPE(TimeGraph);

  private:
    ModelGraph init_;

    ModelGraph main_;

    ModelGraph post_;

    int current_timestep_;

    Scalar current_time_;

    std::string name_;

    Scalar endTime_;

    Scalar deltaT_;

    Scalar max_deltaT_;

    Scalar min_deltaT_;

    int iterations_;

    bool iter_mode = false;

    std::map<std::string, Scalar> model_timestep_restrictions_;

  public:
    TimeGraph(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();

    void execute_pre() { init_.execute(); }

    void execute_main();

    void set_model_timestep(std::string name, Scalar dt) {
        log().info() << "Setting timestep limit " << dt << " for Model "
                     << name;
        model_timestep_restrictions_[name] = dt;
    };

    void execute_post() { post_.execute(); }

    void push_back_pre(std::shared_ptr<Model> m) {
        init_.sub_model_push_back(m);
    }

    void push_back_main(std::shared_ptr<Model> m) {
        main_.sub_model_push_back(m);
    }

    void push_back_post(std::shared_ptr<Model> m) {
        post_.sub_model_push_back(m);
    }

    Scalar get_deltaT();

    float get_maxDeltaT() {
        if (iter_mode) return 1.0;
        return max_deltaT_;
    }

    void set_deltaT(float deltaT) { deltaT_ = deltaT; }

    int &get_current_timestep() { return current_timestep_; }

    void set_iteration(int iteration) { iterations_ = iterations_;}

    void set_current_timestep(int current_timestep) {
        current_timestep_ = current_timestep;
    }
};

#endif
