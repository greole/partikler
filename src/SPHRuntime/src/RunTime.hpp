#ifndef RUNTIME_H
#define RUNTIME_H

#include  <string>
#include  <sstream>
#include  <iostream>
#include <any>

#include "git.hpp"
#include "CGALTYPEDEFS.hpp"
#include "SPHDatastructures.hpp"
#include "FileIO.hpp" // prepare_data_folder

class RunTime {

  private:
    Logger logger_;

    const bool from_disk_;

    // SPHObject& root_object_;
    SPHObjectRegistry obj_reg_;

    size_t n_particles_ = 0;

    float search_cube_size_ = 1.0; // multiple of dx

    // TODO Replace by yaml node
    std::map<std::string, std::any> dict_;

    SearchCubeDomain search_cube_domain_;

    SortedParticles sorted_particles_;

    SortedNeighbours particle_neighbours_;

    Kernel kernel_;

    size_t timestep_ = 0;

    // keep track of particle facet relation
    std::vector<Facet_handle> facets_;
    // std::vector<SPHObject *> fields_;

  public:

    // Constructor

    RunTime(Logger logger, const bool from_disk = true);

    // Setter

    void set_facets(std::vector<Facet_handle> &facets) {
        logger_.info() << " Setting facets ";
        for (auto facet : facets) facets_.push_back(facet);
        logger_.info() << " Done ";
    }

    void set_dict(std::string key, std::any val) { dict_[key] = val; }

    SPHPointField &set_particle_positions(std::vector<Point> &positions);

    RunTime &operator++(int) {
        timestep_++;
        return *this;
    };


    // Getter

    std::vector<Facet_handle> &get_facets() { return facets_; }

    SPHObjectRegistry& get_obj_reg() { return obj_reg_; }

    auto get_dict(std::string key) { return dict_[key]; }

    std::vector<size_t> &get_sorting_idxs() {
        return sorted_particles_.sorting_idxs;
    }

    SPHPointField &get_particle_positions() {
        return  obj_reg_.get_object<SPHPointField &>("Pos");
    }

    size_t get_timestep() const { return timestep_; };

    int get_write_freq() { return std::any_cast<int>(get_dict("write_freq")); }

    bool end() {
        size_t n_timesteps = std::any_cast<int>(dict_["n_timesteps"]);
        return timestep_ >= n_timesteps;
    }

    Logger &log() { return logger_; };

    template <class T> T &get_obj(const std::string name) {
        return obj_reg_.get_object<T &>(name);
    }

    // SPHObject &get_root_object() { return root_object_; }


    // Create

    template <class T> T &create_generic(const std::string name) {
        return obj_reg_.register_object<T>(
            std::make_unique<T>(name, "generic", typename T::value_type()));
    }

    template <class T>
    T &create_generic(const std::string name, typename T::value_type val) {
        return obj_reg_.register_object<T>(
            std::make_unique<T>(name, "generic", val));
    }

    template <class T>
    T &create_field(
        const std::string name,
        typename T::value_type init_value,
        const std::vector<std::string> comp_names) {
        return obj_reg_.register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_, init_value),
            comp_names, name));
    }

    template <class T>
    T &create_field(
        const std::string name, typename T::value_type init_value) {
        return obj_reg_.register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_, init_value),
            name));
    }

    template <class T>
    T &create_field(const std::string name) {
        return obj_reg_.register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_), name));
    }

    SPHSizeTField &create_idx_field() {
        // TODO move to cpp
        auto &f = create_field<SPHSizeTField>("idx");
        std::cout << " n_particles_ " << n_particles_ << std::endl;

        for (size_t i = 0; i < n_particles_; i++) { f[i] = i; }

        return f;
    }

    SortedNeighbours &initialize_particle_neighbours();

    void update_neighbours();

    void update() { n_particles_ = get_particle_positions().size(); }

    Kernel &initialize_kernel();

    void print_timestep() {
        std::stringstream foo;
        foo << get_timestep();
        logger_.info() << "Timestep ";
    };

    // TODO refactor, go through a list of IOObjects and call write method
    void write_to_disk();

    void print_git_state();
};

#endif
