/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

#include "RunTime.hpp"
#include "SearchCubes.cpp"

RunTime::RunTime(Logger logger, const bool from_disk)
    : logger_(logger),
      from_disk_(from_disk),
      obj_reg_({})
{
    logger_.set_scope("Runtime.h");
    logger_.info() << "Starting SPH RunTime ";
    print_git_state();
};

Kernel &RunTime::initialize_kernel() {

    logger_.info_begin() << "Initialising Kernel";

    // const size_t n_pairs = search_cubes_.get_particle_distances().origId.size();

    kernel_ = Kernel {};
    kernel_.W = std::vector<float>(); // (n_pairs, 0);
    kernel_.W.reserve(particle_neighbours_.ids.size());
    kernel_.dWdxo = std::vector<Vector>();
    kernel_.dWdxn = std::vector<Vector>();
    kernel_.dWdxo.reserve(particle_neighbours_.ids.size());
    kernel_.dWdxn.reserve(particle_neighbours_.ids.size());
    logger_.info_end();

    return kernel_;
}

void RunTime::write_to_disk() {
    std::cout << "write to disk" << std::endl;
    int cur_timestep = get_timestep();
    int write_freq = get_write_freq();
    int index_on_dist = cur_timestep / write_freq;
    std::cout << "TIMESTEP: " << cur_timestep << " " << write_freq << std::endl;

    if (cur_timestep - index_on_dist * write_freq == 0) {
        std::cout << obj_reg_.get_objects().size() << std::endl;

        prepare_data_folder("Data", index_on_dist);

        const std::string path = "Data/step#" + intToStr(index_on_dist);

        for (auto &f : obj_reg_.get_objects()) f->write_to_disk(path);
    }
}

searchcubes::SortedNeighbours &
RunTime::initialize_particle_neighbours() {
    // // TODO call set_particle_positions first
    // logger_.info_begin() << "Sorting particles";
    // float dx = std::any_cast<float>(dict_["dx"]);

    // SPHPointField &pos = get_particle_positions();
    // auto & facets  = get_obj<SPHField<Facet_handle>>("facets").get_field();

    // search_cube_domain_ = initSearchCubeDomain(pos.get_field(), search_cube_size_*dx);
    // sorted_particles_ = countingSortParticles(search_cube_domain_, pos.get_field());

    // logger_.info_end();
    // // set_particle_positions(sorted_particles_.particles);

    // reorder_vector(get_sorting_idxs(), facets);

    // logger_.info_begin() << "Initialising particle neighbours";
    // particle_neighbours_ =
    //     createNeighbours(search_cube_domain_, sorted_particles_, facets);
    // logger_.info_end();

    // return particle_neighbours_;
}

void
RunTime::update_neighbours() {
  // // TODO call set_particle_positions first
  // logger_.info_begin() << "Updating particles";
  // float dx = std::any_cast<float>(dict_["dx"]);

  // SPHPointField &pos = get_particle_positions();

  // auto & facets  = get_obj<SPHField<Facet_handle>>("facets").get_field();

  // search_cube_domain_ =
  //     initSearchCubeDomain(pos.get_field(), search_cube_size_ * dx);

  // logger_.info_begin() << "Start counting sort";
  // sorted_particles_ = countingSortParticles(search_cube_domain_, pos.get_field());
  // logger_.info_end();

  // get_particle_positions().set_field(sorted_particles_.particles);

  // logger_.info() << "reoder facets";
  // reorder_vector(get_sorting_idxs(), facets);
  // logger_.info() << "done";

  // logger_.info_begin() << "Updating particle neighbours";
  // particle_neighbours_ =
  //     createNeighbours(search_cube_domain_, sorted_particles_, facets);
  // logger_.info_end();

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
