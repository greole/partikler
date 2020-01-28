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
#include "yaml-cpp/yaml.h"
#include "ObjectRegistry.hpp"   // for ObjectRegistry

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

    float current_time_;

    std::string name_;

    float endTime_;

    float deltaT_;

    float max_deltaT_;

    int iterations_;

    bool iter_mode = false;

  public:
    // TODO Move implementation to cpp file
    TimeGraph(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg)
        : Model(model_name, parameter, objReg),
          init_(ModelGraph("pre", parameter, objReg)),
          main_(ModelGraph("main", parameter, objReg)),
          post_(ModelGraph("post", parameter, objReg)), current_timestep_(0),
          current_time_(0), name_(read_coeff<std::string>("name")),
          endTime_(read_or_default_coeff<float>("endTime", -1.0)),
          deltaT_(read_or_default_coeff<float>("deltaT", -1.0)),
          max_deltaT_(read_or_default_coeff<float>("max_deltaT", -1.0)),
          iterations_(read_or_default_coeff<int>("iterations", 0)) {
        if (deltaT_ < 0) {
            iter_mode = true;
        }
    };

    void execute() {
        execute_pre();
        execute_main();
        execute_post();
    };

    void execute_pre() { init_.execute(); }

    void execute_main() {
        // TODO register some kind of call back
        while (current_timestep_ < iterations_) {
            main_.execute();
            // TODO dont write via objReg, write via Writer SubModel
            // get_objReg().write_to_disk(current_timestep_, name_);
            current_timestep_++;
        };
    }

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

    float get_deltaT() {
        if (iter_mode) return 1.0e-32;
        return deltaT_;
    }

    float get_maxDeltaT() {
        if (iter_mode) return 1.0;
        return max_deltaT_;
    }

    void set_deltaT(float deltaT) { deltaT_ = deltaT; }

    int &get_current_timestep() { return current_timestep_; }
};

#endif
