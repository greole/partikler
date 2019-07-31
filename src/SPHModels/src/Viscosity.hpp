
#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"


class Viscosity : public SPHModel {

    REGISTER_DEC_TYPE(Viscosity);

private:

    // In
    const SPHField<NeighbourPair> &np_;
    const SPHField<VectorPair> &dW_;
    const SPHVectorField &u_;
    const SPHPointField &pos_; // Particle positions

    // Out
    SPHVectorField &dnu_;

    // Coeffs
    float nu_;

public:
    Viscosity(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
