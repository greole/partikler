# include "Viscosity.hpp"


Viscosity::Viscosity(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime),
      // rho_(get_runTime().get_obj<SPHFloatField>("rho")),
      np_(get_runTime().get_obj<SPHField<NeighbourPair>>("neighbour_pairs")),
      dW_(get_runTime().get_obj<SPHField<VectorPair>>("KerneldWdx")),
      u_(get_runTime().get_obj<SPHVectorField>("u")),
      pos_(get_runTime().get_obj<SPHPointField>("Pos")),
      dnu_(get_runTime().get_obj<SPHVectorField>("dnu")),
      nu_(read_or_default_coeff<float>("nu", 1e-05))
{};

void Viscosity::execute() {

    log().info_begin() << "Computing dnu";

    // const SPHFloatField tmp0 = nu.add_ab(pn) / rho.mult_ab(pn);

    const SPHVectorField dxp = particle_distance_vec(pos_, np_);

    // const SPHFloatField tmp1 = (u.sub_ab(pn) * dxp) / dxp.norm();

    // const SPHFloatField tmp = tmp0*tmp1;

    const SPHFloatField tmp = (u_.sub_ab(np_) * dxp)/dxp.norm();

    // TODO Reset 
    dnu_.set(Vector {0,0,0});

    // // weighted sum
    // const size_t size = pn.ids.size();
    // for (size_t i = 0; i < size; i++) {
    //     size_t ownId = pn.ids[i].ownId;
    //     size_t neighId = pn.ids[i].neighId;

    //     for (int j = 0; j < 3; j++) {
    //         dnu[ownId][j]   -= tmp[i]  * kernel.dWdxo[i][j];
    //         dnu[neighId][j] -= tmp[i]  * kernel.dWdxn[i][j];
    //     }
    // }

    dnu_.weighted_sum(np_, tmp, dW_);

    dnu_ *= nu_;

    log().info_end();

};

REGISTER_DEF_TYPE(TRANSPORTEQN, Viscosity);
