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
