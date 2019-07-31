#ifndef PARTICLENEIGHBOURS_H
#define PARTICLENEIGHBOURS_H

#include "SPHModels.hpp"
#include "SPHDatastructures.hpp"

#include "SearchCubes.hpp"

class SPHSTLParticleNeighbours: public SPHModel {

    REGISTER_DEC_TYPE(SPHSTLParticleNeighbours);

private:

    // Coeffs
    float dx_;

    // In
    SPHPointField &pos_;

    SPHField<Facet_handle> & facets_;

    SPHField<SearchCube> & sc_;

    // Out
    SPHField<NeighbourPair> &np_;

    SPHField<STLSurfaceDist> &sd_;

    // Regular data member
    SPHGeneric<SearchCubeDomain> & scd_;

    float search_cube_size_ = 1.0;
public:

    SPHSTLParticleNeighbours(
        const std::string &model_name,
        YAML::Node parameter,
        RunTime & runTime);

    void execute();

    void update_search_cube_domain();

};

#endif
