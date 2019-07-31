
#ifndef MOMENTUM_H
#define MOMENTUM_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"


class Momentum : public SPHModel {

    REGISTER_DEC_TYPE(Momentum);

private:

    // In
    const SPHVectorField &dnu_;
    const SPHVectorField &dp_;

    // Out
    SPHVectorField &u_;
    SPHVectorField &du_;
    SPHGeneric<TimeInfo>& time_;

public:
    Momentum(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();

};

#endif
