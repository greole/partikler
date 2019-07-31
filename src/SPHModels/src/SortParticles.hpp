#ifndef SORTPARTRICLES_H
#define SORTPARTRICLES_H

#include "SPHDatastructures.hpp"
#include "SPHModels.hpp"

#include "SearchCubes.hpp"

class CountingSortParticles : public SPHModel {

    REGISTER_DEC_TYPE(CountingSortParticles);

  private:
    // In
    SPHPointField &pos_;

    SPHField<SearchCube> &sc_;

    //  Sorting indexes
    SPHSizeTField &si_;

    SPHGeneric<SearchCubeDomain> &scd_;

  public:
    CountingSortParticles(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();

    void reorder_fields();
};

#endif
