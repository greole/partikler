#ifndef KERNEL_H
#define KERNEL_H

#include "SPHModels.hpp"
#include "SPHDatastructures.hpp"

class Conti : public SPHModel {

    REGISTER_DEC_TYPE(Conti);

private:
    // In
    const SPHPointField &pos_; // Particle positions

    // const SortedNeighbours &pn_; // Particle neighbours
    const SPHField<NeighbourPair> &np_;
    const SPHFloatField &W_;

    // Out
    SPHFloatField &rho_;

    // Coeffs
    const float lower_limit_;


public:
    Conti(
        const std::string &model_name, YAML::Node parameter, RunTime &runTime);

    void execute();
};

#endif

