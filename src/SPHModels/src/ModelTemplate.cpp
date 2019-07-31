# include "ModelTemplate.hpp"


Model::Model(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime),
{};

void Model::execute() {};

REGISTER_DEF_TYPE(NAMESPACE, Constructor);
