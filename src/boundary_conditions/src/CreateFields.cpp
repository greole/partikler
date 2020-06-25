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

#include "CreateFields.hpp"

#include "Field.hpp" // for PointField
#include "ObjectRegistry.hpp"

InitFields::InitFields(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg) {
    if (parameter["ScalarFields"]) {
        for (auto p : parameter["ScalarFields"])
            float_fields_.push_back(p.as<std::string>());
    }

    if (parameter["VectorFields"]) {
        for (auto p : parameter["VectorFields"])
            vec_fields_.push_back(p.as<std::string>());
    }
}

void InitFields::execute() {
    for (auto f : float_fields_) {
        std::cout << "creating field " << f << std::endl;
        get_objReg().create_field<ScalarField>(f, 0.0);
    }
}

REGISTER_DEF_TYPE(FIELDS, InitFields);
