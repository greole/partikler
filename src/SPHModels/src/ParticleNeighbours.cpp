#include "ParticleNeighbours.hpp"

SPHSTLParticleNeighbours::SPHSTLParticleNeighbours(
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)
    : SPHModel(model_name, parameter, runTime),
      dx_(parameter["dx"].as<float>()),
      pos_(get_runTime().get_particle_positions()),
      facets_(get_runTime().get_obj<SPHField<Facet_handle>>("facets")),
      sc_(get_runTime().get_obj<SPHField<SearchCube>>("search_cubes")),

      np_(get_runTime().get_obj<SPHField<NeighbourPair>>("neighbour_pairs")),
      sd_(get_runTime().get_obj<SPHField<STLSurfaceDist>>("surface_dist")),

      scd_(get_runTime().get_obj<SPHGeneric<SearchCubeDomain>>("search_cube_domain"))
{

    auto pcs = SPHModelFactory::createInstance(
            "SORTING",
            "CountingSortParticles",
            "sort",
            parameter,
            runTime);

    sub_model_push_back(*pcs);
};

void SPHSTLParticleNeighbours::execute() {

    // log().set_scope("")
    log().info() << " Constructing particle neighbours";
    update_search_cube_domain();

    execute_submodels();

    // TODO move sorting to countingSortParticles submodel
    // reorder_vector(
    //     get_runTime().get_obj<SPHSizeTField>("sorting_idxs").get_field(),
    //     facets_.get_field());

    log().info() << "Call createSTLNeighbours";

    // TODO to much copying
    auto [np, sd] = createSTLNeighbours(
        scd_(),
        pos_.get_field(),
        sc_.get_field(),
        facets_.get_field());

    np_.set_field(np);
    sd_.set_field(sd);

    log().info() << "Found ";
    std::cout << np_.size() << " Neighbour Pairs" << std::endl;

    log().info() << "DONE";
    // logger_.info_end();
};

void SPHSTLParticleNeighbours::update_search_cube_domain() {
    scd_() = initSearchCubeDomain(
        pos_.get_field(), search_cube_size_*dx_);
}

REGISTER_DEF_TYPE(PARTICLENEIGHBOURS, SPHSTLParticleNeighbours);
