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

#include <vector>
#include "Models.hpp"

ModelFactory::map_type *ModelFactory::map_;

Equation::Equation(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      np_(objReg.get_object<NeighbourFieldAB>("neighbour_pairs")),
      W_(objReg.get_object<FloatField>("KernelW")),
      dW_(objReg.get_object<KernelGradientField>("KerneldWdx")) {}

// ModelRegister<SPHModelGraph> SPHModelGraph::reg("Core::SPHModelGraph");

REGISTER_DEF_TYPE(CORE, ModelGraph);
REGISTER_DEF_TYPE(CORE, TimeGraph);
