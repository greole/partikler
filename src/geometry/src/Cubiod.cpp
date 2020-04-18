
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

#include "Cubiod.hpp"

#include "Helper.hpp"
#include "Scalar.hpp"
#include "ReaderBase.hpp"

InitShape::InitShape(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ParticleGeneratorBase(model_name, parameter, objReg),
      shape_(read_coeff<std::string>("shape")),
      dimensions_(read_vec3("dimensions")),
      position_(read_vec3("position")),
      noise_(read_or_default_coeff("noise", 0.1)) {}

void InitShape::execute() {

    auto time = get_objReg().get_object<TimeGraph>("TimeGraph");

    if (time.get_current_timestep() == 0) {

        auto cube_pos =
            create_uniform_particle_cube(dimensions_, position_, dx_, noise_);
        // pos_.insert(points_.end(), cube_pos.begin(), cube_pos.end());
        for (auto p : cube_pos)
            pos_.push_back({(Scalar)p[0], (Scalar)p[1], (Scalar)p[2]});

        post_execute();
    }
}

REGISTER_DEF_TYPE(FIELDS, InitShape);
