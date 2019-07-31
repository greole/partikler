#include "Conti.hpp"

Conti::Conti (
    const std::string &model_name, YAML::Node parameter, RunTime &runTime):
    SPHModel(model_name, parameter, runTime),
    pos_(get_runTime().get_particle_positions()),
    np_(get_runTime().get_obj<SPHField<NeighbourPair>>("neighbour_pairs")),
    W_(get_runTime().get_obj<SPHFloatField>("KernelW")),
    rho_(get_runTime().get_obj<SPHFloatField>("rho")),
    lower_limit_(read_or_default_coeff<float>("lower_limit", 0.0))
{};

void Conti::execute() {

    log().info_begin() << "Computing density";

    rho_.set_uniform(0.0);

    rho_.weighted_sum(np_, W_);

    rho_.lower_limit(lower_limit_);

    log().info_end();
};

REGISTER_DEF_TYPE(TRANSPORTEQN, Conti);
