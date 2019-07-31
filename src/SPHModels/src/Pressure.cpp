# include "Pressure.hpp"



Pressure::Pressure(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime),
      rho_(get_runTime().get_obj<SPHFloatField>("rho")),
      np_(get_runTime().get_obj<SPHField<NeighbourPair>>("neighbour_pairs")),
      W_(get_runTime().get_obj<SPHFloatField>("KernelW")),
      dW_(get_runTime().get_obj<SPHField<VectorPair>>("KerneldWdx")),
      p_(get_runTime().get_obj<SPHFloatField>("p")),
      dp_(get_runTime().get_obj<SPHVectorField>("dp")),
      c_(read_or_default_coeff<float>("c", 300.0)),
      rho_0_(read_or_default_coeff<float>("rho_0", 1.0)),
      gamma_(read_or_default_coeff<float>("gamma", 1.4)),
      p_0_(read_or_default_coeff<float>("p_0", 10000)),
      prefac_(c_*c_*rho_0_/gamma_)
{};

void Pressure::execute() {

    log().info_begin() << "Computing pressure";

    const SPHScalarField tmp0 = (rho_/rho_0_).pow(gamma_);
    const SPHScalarField tmp1 = ((tmp0 - 1.0)*prefac_)+p_0_;

    p_ = tmp1;
    log().info_end();

    log().info_begin() << "Computing gradient";

    const SPHFloatField prho = p_/(rho_*rho_);
    const SPHFloatField tmp_ab = prho.add_ab(np_);

    // reset
    dp_.set_uniform({0,0,0});

    dp_.weighted_sum(np_, tmp_ab, dW_);

    log().info_end();
};

REGISTER_DEF_TYPE(TRANSPORTEQN, Pressure);
