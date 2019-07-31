
#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"


class Model : public SPHModel {

    REGISTER_DEC_TYPE(Model);

private:


public:
    Model(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
