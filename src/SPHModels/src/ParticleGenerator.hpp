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

#include "SPHModels.hpp"
#include "SPHDatastructures.hpp"

#include <memory>
// TODO needs a polyhedron generator

class SPHSTLReader: public SPHModel {

    // Reads an stl file and creates a Polyhedron representation
    REGISTER_DEC_TYPE(SPHSTLReader);

private:

    std::string fn_;


public:

   SPHSTLReader(
        const std::string &model_name,
        YAML::Node parameter,
        RunTime & runTime);

    void execute();

};


class SPHParticleGenerator: public SPHModel {


    // A Polyhedron particle generator
    //
    // Dependencies CGALPolyhedron, facets, pos
    //
    REGISTER_DEC_TYPE(SPHParticleGenerator);

private:

    // The base polyhedron on which particles are created
    SPHGeneric<CGALPolyhedron>& polyhedron_;

    SPHField<Facet_handle> & facets_;

    // Reference to positions field
    SPHPointField& pos_;

    const float dx_;


public:

    SPHParticleGenerator(
        const std::string &model_name,
        YAML::Node parameter,
        RunTime & runTime);


    void execute();

};

#endif
