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


#ifndef PARTIKLER_CREATEFIELDS_INCLUDED_H
#define PARTIKLER_CREATEFIELDS_INCLUDED_H

#include "Models.hpp"
#include "cgal/CGALHelper.hpp"

class InitFields: public Model {

    REGISTER_DEC_TYPE(InitFields);

private:

    PointField& pos_;

    std::vector<std::string> float_fields_ {};

    std::vector<std::string> vec_fields_ {};


public:

InitFields(
    const std::string &model_name,
    YAML::Node parameter,
    ObjectRegistry & objReg);

    void execute();

};

#endif
