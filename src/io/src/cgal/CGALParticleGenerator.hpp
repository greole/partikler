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

#ifndef PARTICLE_GENERATOR_H
#define PARTICLE_GENERATOR_H

#include <fstream>
#include <memory>
#include <string> // for string
#include <vector> // for vector

#include "Field.hpp"  // for IntField, Field (ptr only), PointField
#include "Models.hpp" // for Model, ModelRegister (ptr only)
#include "cgal/CGALHelper.hpp"
#include "cgal/CGALTYPEDEFS.hpp" // for CGALPolyhedron, Facet_handle

class ObjectRegistry;

namespace YAML {
class Node;
} // namespace YAML
template <class T> class Generic;
// TODO needs a polyhedron generator

class SPHSTLReader : public Model {

    // Reads an stl file and creates a Polyhedron representation
    REGISTER_DEC_TYPE(SPHSTLReader);

  private:
    std::string fn_;

  public:
    SPHSTLReader(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

class SPHParticleGenerator : public Model {

    // A Polyhedron particle generator
    //
    // Dependencies CGALPolyhedron, facets, pos
    //
    REGISTER_DEC_TYPE(SPHParticleGenerator);

  private:
    // FieldIdMap& fieldIdMap_;

    // id of the boundary
    int boundary_id_;

    // The base polyhedron on which particles are created
    Generic<CGALPolyhedron> &polyhedron_;

    Field<std::vector<Facet_handle>> &facets_;

    VectorField &pos_;

    SizeTField &idx_;

    IntField &type_;

    IntField &boundary_;

    const float dx_;

  public:
    SPHParticleGenerator(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
