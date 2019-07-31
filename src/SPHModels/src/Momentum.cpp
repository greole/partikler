# include "Momentum.hpp"


Momentum::Momentum(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime),
      dnu_(get_runTime().get_obj<SPHVectorField>("dnu")),
      dp_(get_runTime().get_obj<SPHVectorField>("dp")),
      u_(get_runTime().get_obj<SPHVectorField>("u")),
      du_(get_runTime().get_obj<SPHVectorField>("du")),
      time_(get_runTime().get_obj<SPHGeneric<TimeInfo>>("TimeInfo"))
{};

void Momentum::execute() {

    log().info_begin() << "Computing du/dt";

    du_ = dnu_ - dp_;

    log().info_end();

    log().info_begin() << "Computing velocity";

    // TODO implement time integrator
     u_ = u_/1.3;

    // CFL = max(u)*deltaT/dx

    float maxU = u_.norm().get_max();

    std::cout << "maxU"
              << maxU
              << "deltaT" <<  time_().deltaT
              << "maxCFL"
              << maxU * time_().deltaT /
                     std::any_cast<float>(get_runTime().get_dict("dx"))
              << std::endl;

    // std::cout << du_ << std::endl;
    // std::cout << du_*time_().deltaT << std::endl;
    // std::cout << u_ << std::endl;

    u_ += (du_ * time_().deltaT);

    log().info_end();
};

REGISTER_DEF_TYPE(TRANSPORTEQN, Momentum);
