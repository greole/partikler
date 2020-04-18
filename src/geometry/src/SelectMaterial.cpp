
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

#include "SelectMaterial.hpp"

#include "Helper.hpp"
#include "Scalar.hpp"
#include "ReaderBase.hpp"

SelectMaterial::SelectMaterial(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      fieldIdMap_(objReg.get_object<FieldIdMap>("FieldIdMap")),
      materialMap_(objReg.get_object<MaterialMap>("MaterialMap")),
      material_(materialMap_.getMaterial(read_coeff<std::string>("material"))),
      name_(read_coeff<std::string>("name")),
      fieldId_(fieldIdMap_.append(name_, material_)),
      pos_(objReg.create_field<VectorField>("Pos", {}, {"X", "Y", "Z"})),
      id_(objReg.create_field<IntField>("id")),
      shape_(read_coeff<std::string>("shape")),
      start_(read_vec3("start")),
      end_(read_vec3("end"))
{}

// check if pos is inside box
bool is_inside_box(Vec3 &pos, std::pair<Vec3, Vec3> box) {
    if (pos[0] > box.first[0] && pos[1] > box.first[1] &&
        pos[2] > box.first[2] && pos[0] < box.second[0] &&
        pos[1] < box.second[1] && pos[2] < box.second[2]) {
        return true;
    }
    return false;
}

void SelectMaterial::execute() {


    // get existing particles
    // check if inside box
    std::pair<Vec3, Vec3> bb {start_,end_};
    for(size_t i=0; i<pos_.size(); i++) {
        if ( is_inside_box(pos_[i], bb)) {
                // change id and material
            id_[i] = fieldId_;
        }

    }

}

REGISTER_DEF_TYPE(FIELDS, SelectMaterial);
