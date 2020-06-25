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

#include "ForcesList.hpp"

#include "Time.hpp"

ForcesList::ForcesList(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "fBody", {}, {"fBodyx", "fBodyy", "fBodyz"})) {

    for (auto el : parameter["submodel"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();
        auto params = it->second;

        auto model = ModelFactory::createInstance(
            model_namespace, model_name, model_name, params, objReg);

        sub_model_push_back(model);
    }
}

void ForcesList::execute() {

    log().info_begin() << "Collecting Forces";

    // reset
    for (size_t i = 0; i < f_.size(); i++) {
        f_[i] = {0.0, 0.0, 0.0};
    }

    for (auto &force : get_submodels()) {
        std::cout << force.get()->get_name() << std::endl;
        dynamic_cast<VectorFieldEquation *>(force.get())->execute();
        auto &fm = dynamic_cast<VectorFieldEquation *>(force.get())->get();
        for (size_t i = 0; i < f_.size(); i++) {
            f_[i] += fm[i];
        }
    }

    log().info_end();
}

REGISTER_DEF_TYPE(FORCES, ForcesList);
