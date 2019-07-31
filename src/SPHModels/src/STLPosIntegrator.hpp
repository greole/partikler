
#ifndef STLPOSINTEGRATOR_H
#define STLPOSINTEGRATOR_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"
#include "STLLimitedDx.hpp"

class STLPosIntegrator : public SPHModel {

    REGISTER_DEC_TYPE(STLPosIntegrator);

private:

    // const std::vector<Point> opoints_;
    const SPHIntField &type_;
    SPHField<Facet_handle> & facets_;
    const SPHSizeTField &idx_;

    // Out
    SPHVectorField &u_;
    SPHPointField& pos_;
    SPHGeneric<TimeInfo>& time_;

public:
    STLPosIntegrator(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
