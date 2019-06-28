#ifndef RUNTIME_H
#define RUNTIME_H

#include  <string>
#include  <sstream>
#include  <iostream>

#include "git.hpp"
#include "SPHio.hpp"
#include <any>

class RunTime {

  private:
    Logger logger_;

    const bool from_disk_;

    SPHObject root_object_;

    bool constructed_from_disk_;

    size_t n_particles_;

    float search_cube_size_ = 1.0; // multiple of dx

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

    RunTime(Logger logger, const bool from_disk = true)
      : logger_(logger), from_disk_(from_disk), root_object_("root", "SPHObject") {
        logger_.set_scope("Runtime.h");
        logger_.info() << "Starting SPH RunTime ";
        print_git_state();

        if (from_disk_) {
            std::cout << "[INFO] Reading input files from disk " << std::endl;

            constructed_from_disk_ = true;
        }
    };

    void set_facets(
        std::vector<Facet_handle>& facets
        ) {
        logger_.info() << " Setting facets ";
        for (auto facet: facets)  facets_.push_back(facet);
        logger_.info() << " Done ";
    }

    std::vector<Facet_handle>& get_facets(
        ) {return facets_;}


    void set_dict(std::string key, std::any val) { dict_[key] = val; }

    auto get_dict(std::string key) {return dict_[key]; }

  std::vector<size_t> & get_sorting_idxs() {return sorted_particles_.sorting_idxs;}


    // TODO use size_t here
    std::vector<int> get_searchCubeIds() const {
        std::vector<int> ret(n_particles_, 0);

        size_t sc_ctr = 0;
        // for (auto search_cube: search_cubes_.get_searchCubes()) {
        //     for(size_t p_ctr = search_cube.first; p_ctr <= search_cube.last;
        //     p_ctr++) {
        //         ret[p_ctr] = (int) sc_ctr;
        //     }
        //     sc_ctr++;
        // }

        return ret;
    }

    // SPHIntField& create_uniform_int_field(
    //         const int init_value,
    //         const std::string name
    //         ) {

    //     // TODO call destructors or use some sort of smart_ptr
    //     std::vector<int> * f = new std::vector<int> (n_particles_,
    //     init_value); SPHIntField* SPHField = new SPHIntField(*f, name);
    //     register_field(*SPHField);

    //     return *SPHField;
    // }

    SPHIntField& create_uniform_int_field (
        const int init_value,
        const std::string name
        ) {
        // TODO call destructors or use some sort of smart_ptr
        logger_.info() << "Creating int field " << name;
        std::vector<int> *f =
            new std::vector<int>(n_particles_, init_value);
        SPHIntField *SPHField = new SPHIntField(*f, name);
        // root_object_.register_object(*SPHField);

        return *SPHField;
    }


    SPHFloatField& create_uniform_field (
            const float init_value,
            const std::string name
            ) {
        // TODO call destructors or use some sort of smart_ptr
        logger_.info() << "Creating float field " << name;
        std::vector<float> *f =
            new std::vector<float>(n_particles_, init_value);
        SPHFloatField *SPHField = new SPHFloatField(*f, name);
        // root_object_.register_object(*SPHField);

        return *SPHField;
    }

    SPHVectorField &create_uniform_field(
        const Vector init_value,
        const std::string name,
        const std::vector<std::string> comp_names) {
        logger_.info() << "Creating vector field " << name;
        std::vector<Vector> *f =
            new std::vector<Vector>(n_particles_, init_value);

        SPHVectorField *SPHField = new SPHVectorField(*f, comp_names, name);
        // root_object_.register_object(*SPHField);

        return *SPHField;
    }

  SPHSizeTField& create_idx_field(){
    logger_.info() << "Creating idx field ";
    std::vector<size_t> *f =
      new std::vector<size_t>(n_particles_, 0);

    for(size_t i=0;i<n_particles_;i++){
      f->operator[](i)=i;
    }

    SPHSizeTField *SPHField = new SPHSizeTField(*f, "idx");
    // root_object_.register_object(*SPHField);

    return *SPHField;
  }

    SPHPointField &set_particle_positions(std::vector<Point> &positions);

    SortedNeighbours &
    initialize_particle_neighbours(std::vector<Point> &positions);

    void update_neighbours();

    SPHPointField& get_particle_positions() {
        return root_object_.get_object<SPHPointField&>("Pos");
    }

    Kernel& initialize_kernel();

    size_t get_timestep() const {return timestep_;};

    RunTime& operator++(int) {
        timestep_++;
        return *this;
    };

    void print_timestep() {
      std::stringstream foo;
      foo <<  get_timestep();
      logger_.info() << "Timestep ";
    };

    // ParticleNeighbours& update_particle_neighbours(std::vector<Point>&
    // sorted_points) {
    //     return search_cubes_.update_particle_neighbours(sorted_points);
    // };

    // SearchCubeTree & get_search_cubes() {return search_cubes_;};

    // size_t get_timeStep() const {return timestep_;};

    // void register_field(SPHObject &f) {
    //     std::cout << "[INFO] Registering field: " << f.get_name() << std::endl;
    //     fields_.push_back(&f);
    // }

    // template< class T>
    // T& get_field(const std::string name) const {
    //     for (SPHObject* f: fields_) {
    //         if (f->get_name() == name) {return dynamic_cast<T&> (*f);};
    //     }
    // }

  int get_write_freq(){ 
    return std::any_cast<int>(get_dict("write_freq"));
  }

  // TODO refactor, go through a list of IOObjects and call write method
  void write_to_disk() {
    int cur_timestep = get_timestep();
    int write_freq = get_write_freq();
    int index_on_dist = cur_timestep/write_freq;
    std::cout << "TIMESTEP: " <<  cur_timestep << " " << write_freq << std::endl;

    if (cur_timestep - index_on_dist*write_freq == 0) {
       std::cout << root_object_.get_objects().size() << std::endl;

        for (const SPHObject *f : root_object_.get_objects()) {
            std::string data_type = f->get_type();
            std::cout << data_type << " " << f->get_name() << std::endl;

            auto write_field_wrapper = [&](auto &field, std::string length="32") {
                write_field(
                    "Data",
                    index_on_dist,
                    field.get_name(),
                    field.get_field(),
                    field.get_type(),
                    length);
            };

            if (data_type == "int") {
                const SPHIntField &field = static_cast<const SPHIntField &>(*f);
                write_field_wrapper(field);
            };

            if (data_type == "long") {
              // const SPHSizeTField &field = static_cast<const SPHSizeTField &>(*f);
              // // convert to int because of SuperSPHReader
              // std::vector<float> vec (field.get_field().size(), 0);
              // for(size_t i=0;i<field.get_field().size();i++){
              //   vec[i] = (float)field.get_field()[i];
              // }

              // SPHFloatField conv(vec, "idx");
              // write_field_wrapper(conv);
            };

            if (data_type == "float") {
                const SPHFloatField &field =
                    static_cast<const SPHFloatField &>(*f);
                write_field_wrapper(field);
            };

            if (data_type == "vector") {
                const SPHVectorField &field = static_cast<const SPHVectorField &>(*f);

                write_vector_field(
                    "Data",
                    index_on_dist,
                    f->get_name(),
                    field.get_field(),
                    field.get_type(),
                    "32",
                    field.get_comp_names());
            };

            if (f->get_type() == "Point") {
              const SPHPointField &field = static_cast<const SPHPointField &>(*f);

              write_point_field(
                                 "Data",
                                 index_on_dist,
                                 f->get_name(),
                                 field.get_field(),
                                 field.get_type(),
                                 "32",
                                 field.get_comp_names());
            };
        }
    }

    };

    bool end() {
        size_t n_timesteps = std::any_cast<int>(dict_["n_timesteps"]);
        return timestep_ >= n_timesteps;
    }

    void print_git_state();
};

#endif
