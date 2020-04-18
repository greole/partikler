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

#include "ObjectRegistry.hpp"

bool ObjectRegistry::object_exists(const std::string name) const {
    auto search = objects_.find(name);
    if (search != objects_.end()) {
        return true;
    }
    return false;
}

int FieldIdMap::getId(const std::string name) {
    for (size_t i = 0; i < fields_.size(); i++) {
        if (name == fields_[i]) return i;
    }
    std::string error_str = "no field " + name + " found in fieldIdMap";
    throw std::runtime_error(error_str);
}


int FieldIdMap::append(std::string field_name, Material m) {
    // TODO check if name exists and return existing id
    auto res = std::find(fields_.begin(), fields_.end(), field_name);
    if (res != fields_.end()) {
        int id = std::distance(fields_.begin(), res);
        return id;
    }

    // otherwise append new field and material
    int id = fields_.size();
    fields_.push_back(field_name);
    material_.push_back(m);
    return id;
}
