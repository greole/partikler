#include "RunTime.hpp"
#include "SearchCubes.cpp"

Kernel &RunTime::initialize_kernel() {

    logger_.info_begin() << "Initialising Kernel";

    // const size_t n_pairs = search_cubes_.get_particle_distances().origId.size();

    kernel_ = Kernel {};
    kernel_.W = std::vector<float>(); // (n_pairs, 0);
    kernel_.W.reserve(particle_neighbours_.ownId.size());
    kernel_.dWdx = std::vector<Vector>();
    kernel_.dWdx.reserve(particle_neighbours_.ownId.size());
    logger_.info_end();

    return kernel_;
}

SPHPointField &RunTime::set_particle_positions(std::vector<Point> &positions) {
    // creates the positions field from a vector of positions
    // returns a reference to the newly created field
    // TODO delete all registered fields when runTime destructor is called

    SPHPointField *pointField =
      new SPHPointField(sorted_particles_.particles, "Pos");
    register_field(*pointField);

    n_particles_ = sorted_particles_.particles.size();

    return *pointField;
}

SortedNeighbours &
RunTime::initialize_particle_neighbours(std::vector<Point> &positions) {
    // TODO call set_particle_positions first
    logger_.info_begin() << "Sorting particles";
    float dx = std::any_cast<float>(dict_["dx"]);
    search_cube_domain_ = initSearchCubeDomain(positions, 2.8*dx);
    sorted_particles_ = countingSortParticles(search_cube_domain_, positions);

    logger_.info_end();
    set_particle_positions(sorted_particles_.particles);

    logger_.info_begin() << "Initialising particle neighbours";
    particle_neighbours_ =
        createNeighbours(search_cube_domain_, sorted_particles_);
    logger_.info_end();

    return particle_neighbours_;
}

void
RunTime::update_neighbours() {
  // TODO call set_particle_positions first
  logger_.info_begin() << "Updating particles";
  float dx = std::any_cast<float>(dict_["dx"]);

  SPHPointField &pos = get_particle_positions();

  search_cube_domain_ = initSearchCubeDomain(pos.get_field(), 2.8 * dx);

  logger_.info() << "Start counting sort";
  sorted_particles_ = countingSortParticles(search_cube_domain_, pos.get_field());
  logger_.info() << "end counting sort";

  logger_.info_end();

  get_particle_positions().set_field(sorted_particles_.particles);

  logger_.info_begin() << "Updating particle neighbours";
  particle_neighbours_ =
    createNeighbours(search_cube_domain_, sorted_particles_);
  logger_.info_end();

}

void RunTime::print_git_state() {
    // TODO use runtimes info/warn functions
    if (GIT_RETRIEVED_STATE) {
        std::cout << "[INFO] Git commit " << GIT_HEAD_SHA1 << std::endl;
        if (GIT_IS_DIRTY) {
            std::cerr << "[WARN] There were uncommitted changes." << std::endl;
        }
    } else {
        std::cerr << "[WARN] Failed to get the current git state."
                  << "Is this a git repo?" << std::endl;
    }
};
