
#ifndef Pressure_H
#define Pressure_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"

class Pressure : public SPHModel {

    REGISTER_DEC_TYPE(Pressure);

  private:

    // In
    // Density
    const SPHFloatField &rho_;
    const SPHField<NeighbourPair> &np_;
    const SPHFloatField &W_;
    const SPHField<VectorPair> &dW_;

    // Out
    // Pressure
    SPHFloatField &p_;
    SPHVectorField &dp_;

    // Coeffs
    const float c_;
    const float rho_0_;
    const float gamma_;
    const float p_0_;
    const float prefac_;

  public:

    Pressure(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
