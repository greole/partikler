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

#include "Datastructures.hpp"
#include "Models.hpp"
#include "yaml-cpp/yaml.h"

class GenerateBoundaryParticles : public Model {

    REGISTER_DEC_TYPE(GenerateBoundaryParticles);

  private:

    ObjectRegistry local_objReg_;
    TimeGraph& timeGraph_;

    IntField &boundaryIds_;
    IntField &typeIds_;
    SizeTField &idx_;
    PointField &pos_;

    int iterations_;
    int write_freq_;

    std::string filename_;
    std::string boundary_name_;

    float dx_;



  public:

    YAML::Node default_graph();

    GenerateBoundaryParticles(
        const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg);

    void execute();
};

#endif