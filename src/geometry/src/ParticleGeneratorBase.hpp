
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

#ifndef PARTIKLER_PARTICLEGENERATORBASE_INCLUDED_H
#define PARTIKLER_PARTICLEGENERATORBASE_INCLUDED_H

#include "Field.hpp"
#include "Models.hpp"         // for ModelRegister (ptr only), REGISTER_DEC_TYPE
#include "ObjectRegistry.hpp" // for ObjectRegistry
#include "yaml-cpp/yaml.h"

class ParticleGeneratorBase : public Model {

  protected:
    FieldIdMap &fieldIdMap_;

    MaterialMap &materialMap_;

    Material material_;

    ObjectRegistry local_objReg_;

    // PointField &points_;

    VectorField &pos_;

    IntField &id_;

    std::string name_;

    int fieldId_;

    float dx_;

    Vec3 translation_vector_;

  public:
    ParticleGeneratorBase(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    // Vec3 read_vector(YAML::Node parameter, std::string coeff);

    template <class T> void append(T &, std::string name) {
        auto &oreg = get_objReg();
        field_append(
            oreg.get_object<T>(name), local_objReg_.get_object<T>(name));
    }

    virtual void execute() {};

    void post_execute() {

        logger_.info_begin() << "Translate: " << translation_vector_;

        translatePoints(pos_, translation_vector_);

        logger_.info_end();

        // Transfer particles to main object registry
        // Here only the position and the boundary id and type id
        // need to be transferred
        logger_.info_begin() << "Transfering ";

        auto &oreg = get_objReg();
        if (oreg.object_exists("Pos")) {
            append(pos_, "Pos");
        } else {
            // Move the object if it doesn't exist in the main registry yet
            oreg.get_objects()["Pos"] = local_objReg_.get_object_ptr("Pos");
        }

        id_.reserve(pos_.size());
        for (size_t i = 0; i < pos_.size(); i++) {
            id_.push_back(fieldId_);
        }

        if (oreg.object_exists("id")) {
            append(id_, "id");
        } else {
            // Move the object if it doesn't exist in the main registry yet
            oreg.get_objects()["id"] = local_objReg_.get_object_ptr("id");
        }

        get_objReg().update_n_particles();
        logger_.info_end();
    };
};

#endif
