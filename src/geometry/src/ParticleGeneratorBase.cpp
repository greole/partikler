
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

#include "ParticleGeneratorBase.hpp"

ParticleGeneratorBase::ParticleGeneratorBase(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      fieldIdMap_(objReg.get_object<FieldIdMap>("FieldIdMap")),
      materialMap_(objReg.get_object<MaterialMap>("MaterialMap")),
      material_(materialMap_.getMaterial(read_coeff<std::string>("material"))),
      local_objReg_(ObjectRegistry()),
      // points_(local_objReg_.create_field<PointField>(
      //     "Points", {}, {"X", "Y", "Z"})),
      pos_(local_objReg_.create_field<VectorField>("Pos", {}, {"X", "Y", "Z"})),
      id_(local_objReg_.create_field<IntField>("id")),
      name_(read_coeff<std::string>("name")),
      fieldId_(fieldIdMap_.append(name_, material_)),
      dx_(read_coeff<float>("dx")),
      translation_vector_(read_vector(parameter, "translate")) {}

Vec3 ParticleGeneratorBase::read_vector(
    YAML::Node parameter, std::string coeff) {

    if (!parameter[coeff]) return {0, 0, 0};

    auto p = parameter[coeff];
    return {p[0].as<float>(), p[1].as<float>(), p[2].as<float>()};
}
