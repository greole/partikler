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

// TODO: Add .cpp file to CMakeLists.txt

#ifndef GENERATEBOUNDARYPARTICLES_H
#define GENERATEBOUNDARYPARTICLES_H

#include <string> // for string
#include <vector> // for vector
#include <memory> // for vector

#include "Field.hpp" // for PointField
#include "FieldOps.hpp"
#include "Models.hpp"           // for Model, ModelRegister (ptr only), REG...
#include "ObjectRegistry.hpp"   // for FieldIdMap (ptr only), ObjectRegistry
#include "yaml-cpp/node/node.h" // for Node
#include "yaml-cpp/yaml.h"

class GenerateBoundaryParticles : public Model {

    REGISTER_DEC_TYPE(GenerateBoundaryParticles);

  private:
    FieldIdMap &fieldIdMap_;
    ObjectRegistry local_objReg_;
    TimeGraph &timeGraph_;

    // PointField &pos_;

    int iterations_;

    int write_freq_;

    std::string filename_;

    std::string boundary_name_;

    int fieldId_;

    float scale_;

    float dx_;

    Vec3 translation_vector_;

  public:
    YAML::Node default_graph();

    GenerateBoundaryParticles(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    Vec3 read_translation_vector(YAML::Node parameter);

    void execute();

    template <class T> void append(T &, std::string name);
};

#endif
