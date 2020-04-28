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

#include "Equation.hpp"

template<>
template<>
Vec3 FieldEquationBase<VectorField>::read_boundary(YAML::const_iterator it){
    std::cout << "read vector val" << std::endl;
     Vec3 ret {it->second["value"][0].as<Scalar>(),
            it->second["value"][1].as<Scalar>(),
            it->second["value"][2].as<Scalar>()};
     std::cout << ret << std::endl;
     return ret;
}

template<>
template<>
Scalar FieldEquationBase<ScalarField>::read_boundary(YAML::const_iterator it){
    return it->second["value"].as<Scalar>();
}
