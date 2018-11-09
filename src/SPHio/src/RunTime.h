#ifndef RUNTIME_H
#define RUNTIME_H

#include "git.h"
#include <any>

class RunTime {

    private:

        Logger logger_;

        const bool from_disk_;

        bool constructed_from_disk_;

        size_t n_particles_;

        std::map<std::string, std::any> dict_;

        SearchCubeTree search_cubes_;

        // ParticleNeighbours& particle_neighbours_;

        Kernel kernel_;

        size_t timestep_ = 0;

        std::vector<SPHFieldBase*> fields_;

    public:

        // Constructor

        RunTime(
                Logger logger,
                const bool from_disk=true
                ):
            logger_(logger),
            from_disk_(from_disk)
        {
            logger_.set_scope("Runtime.h");
            logger_.info() << "Starting SPH RunTime ";
            print_git_state();

            if (from_disk_){
                std::cout << "[INFO] Reading input files from disk " << std::endl;

                constructed_from_disk_ = true;
            }

        };

        void set_dict(std::string key, std::any val) {

            dict_[key] = val;

        }

        auto get_dict(std::string key) {

            return dict_[key];

        }

        // TODO use size_t here
        std::vector<int>  get_searchCubeIds() const  {
            std::vector<int> ret (n_particles_, 0);

            size_t sc_ctr = 0;
            for (auto search_cube: search_cubes_.get_searchCubes()) {
                for(size_t p_ctr = search_cube.first; p_ctr <= search_cube.last; p_ctr++) {
                    ret[p_ctr] = (int) sc_ctr;
                }
                sc_ctr++;
            }

            return ret;
        }

        SPHIntField& create_uniform_int_field(
                const int init_value,
                const std::string name
                ) {

            // TODO call destructors or use some sort of smart_ptr
            std::vector<int> * f = new std::vector<int> (n_particles_, init_value);
            SPHIntField* SPHField = new SPHIntField(*f, name);
            register_field(*SPHField);

            return *SPHField;
        }


        SPHScalarField& create_uniform_field (
                const float init_value,
                const std::string name
                ) {
            // TODO call destructors or use some sort of smart_ptr
            std::vector<float> * f = new std::vector<float> (n_particles_, init_value);
            SPHScalarField* SPHField = new SPHScalarField(*f, name);
            register_field(*SPHField);

            return *SPHField;
        }

        SPHVectorField& create_uniform_field(
                const Vector init_value,
                const std::string name
                ) {
            std::vector<Vector> * f = new std::vector<Vector> (n_particles_, init_value);
            SPHVectorField* SPHField = new SPHVectorField(*f, name);
            register_field(*SPHField);

            return *SPHField;
        }

        SPHPointField& set_particle_positions(std::vector<Point>& positions) {

            n_particles_ = positions.size();
            SPHPointField* pointField = new SPHPointField(positions, "Pos");
            register_field(*pointField);

            return *pointField;

        }

        SPHPointField& get_particle_positions() {
            return get_field<SPHPointField&>("Pos");
        }

        ParticleNeighbours& initialize_particle_neighbours() {

             float dx = std::any_cast<float>(dict_["dx"]);

             std::vector<Point> & points = get_particle_positions().get_field();

             search_cubes_ = SearchCubeTree (
                     get_particle_positions().get_field(),
                     3*dx, logger_);

             points = search_cubes_.get_sorted_points();

             return search_cubes_.get_particle_distances();
        }

        Kernel& initialize_kernel() {

                logger_.info_begin() << "Initialising Kernel";

                const size_t n_pairs = search_cubes_
                    .get_particle_distances()
                    .origId
                    .size();

                kernel_ = Kernel {};
                kernel_.W    = std::vector<float>  (n_pairs, 0);
                kernel_.dWdx = std::vector<Vector> (n_pairs, {0, 0, 0});
                logger_.info_end();

                return kernel_;
        }

        RunTime& operator++(int) {
            timestep_++;
            return *this;
        };

        ParticleNeighbours& update_particle_neighbours(std::vector<Point>& sorted_points) {
            return search_cubes_.update_particle_neighbours(sorted_points);
        };

        SearchCubeTree & get_search_cubes() {return search_cubes_;};

        size_t get_timeStep() const {return timestep_;};

        void register_field(SPHFieldBase& f) {
            std::cout << "[INFO] Registering field: " << f.get_name() << std::endl;
            fields_.push_back(&f);
        }

        template< class T>
        T& get_field(const std::string name) const {
            for (SPHFieldBase* f: fields_) {
                if (f->get_name() == name) {return dynamic_cast<T&> (*f);};
            }
        }

        void write_to_disk() const {
            for (SPHFieldBase* f: fields_) {
                if (f->get_type() == "Int") {
                    const SPHIntField& field = static_cast<SPHIntField&> (*f);
                    write_field(
                       "Data",
                       timestep_,
                       f->get_name(),
                       field.get_field());
                };
                if (f->get_type() == "Point") {
                    const SPHPointField& field = static_cast<SPHPointField&> (*f);
                    write_field(
                       "Data",
                       timestep_,
                       f->get_name(),
                       field.get_field());
                };
                if (f->get_type() == "Vector") {
                    const SPHVectorField& field = static_cast<SPHVectorField&> (*f);
                    write_field(
                       "Data",
                       timestep_,
                       f->get_name(),
                       field.get_field());
                };
                if (f->get_type() == "Scalar") {
                    const SPHScalarField& field = static_cast<SPHScalarField&> (*f);
                    write_field(
                       "Data",
                       timestep_,
                       f->get_name(),
                       field.get_field());
                };
            }
        }

        bool end() {
            size_t n_timesteps = std::any_cast<size_t>(dict_["n_timesteps"]);
            return timestep_ >= n_timesteps;
        }

        void print_git_state() {
            if(GIT_RETRIEVED_STATE) {
                std::cout
                    << "[INFO] Git commit "
                    << GIT_HEAD_SHA1
                    << std::endl;
                if(GIT_IS_DIRTY) {
                    std::cerr
                        << "[WARN] There were uncommitted changes."
                        << std::endl;
                }
            }
            else {
                std::cerr
                    << "[WARN] Failed to get the current git state."
                    << "Is this a git repo?"
                    << std::endl;
            }
        };



};

#endif
