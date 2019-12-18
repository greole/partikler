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


#ifndef MOMENTUM_H
#define MOMENTUM_H

#include "Models.hpp"
#include "Field.hpp"
#include "FieldOps.hpp"
#include "Vec3.hpp"

#include "yaml-cpp/yaml.h"


class Momentum : public Equation {

    REGISTER_DEC_TYPE(Momentum);

private:

    // In
    const VectorField &dnu_;
    const VectorField &dp_;

    // Out
    VectorField &u_;
    VectorField &du_;

    TimeGraph& time_;

public:
    Momentum(
        const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg);

    void execute();

};

#endif
