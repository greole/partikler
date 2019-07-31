#ifndef KERNEL_H
#define KERNEL_H

#include "SPHModels.hpp"
#include "yaml-cpp/yaml.h"
#include "SPHDatastructures.hpp"

class STLWendland2D : public SPHModel {

    REGISTER_DEC_TYPE(STLWendland2D);

  private:
    // Coeffs
    const float h_;  // Smoothing length
    const float ih_; // Inverse Smoothing length
    // // 3d
    // const float W_fak2 = 21. / (256. * M_PI * h * h * h);
    // const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);
    // 2d
    const float W_fak2_;  // = 7. / (64. * M_PI * h * h);
    const float dW_fak2_; // = 7. / (64. * M_PI * h * h * h);

    // In
    const SPHPointField &pos_; // Particle positions

    const SPHField<NeighbourPair> &np_;
    const SPHField<STLSurfaceDist> &sd_;

    // Out
    // Kernel &kernel                               // Kernel field
    SPHFloatField &W_;
    SPHField<VectorPair> &dWdx_;

  public:
    STLWendland2D(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif
