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

#include "ReaderBase.hpp"


ReaderBase::ReaderBase(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      read_(false),
      time_graph_(objReg.get_object<TimeGraph>("TimeGraph")) {}

HDF5Reader::HDF5Reader(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : ReaderBase("Reader", parameter, objReg),
      file_name_(parameter["name"].as<std::string>()),
      write_freq_(parameter["freq"].as<int>())
{}

void HDF5Reader::execute() {

    const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

    int comm_rank = 0;
    std::string fname {(file_name_ + ".h5part")};

    H5AbortOnError();
    H5SetVerbosityLevel(h5_verbosity);

    // open file and create first step
    h5_prop_t prop = H5CreateFileProp();
    H5SetPropFileCoreVFD(prop, 0);
    h5_file_t file = H5OpenFile(fname.c_str(), H5_O_RDONLY, prop);
    h5_int64_t start, end;

    h5_int64_t step_num = H5GetNumSteps (file);

    std::cout << " Reading file: " << file_name_ << " Timestep: " << step_num  << std::endl;
    H5SetStep(file, 0);
    h5_int64_t n_datasets = H5PartGetNumDatasets (file);
    std::cout << " Number  of datasets " << n_datasets << std::endl;
    h5_int64_t num_particles = H5PartGetNumParticles (file);

    std::cout << " Reading file: " << file_name_ << " Timestep: " << step_num <<
        " Number of particles: " << num_particles << std::endl;

    char name[H5_MAX_NAME_LEN];
    h5_int64_t type;
    h5_size_t dim;

    // NOTE for now read only Position, id, and idx
    // TODO implement a better system to detect fields
    // which are needed for time derivatives

    // Start with position
    for (h5_int64_t i = 0; i < n_datasets; i++) {
        H5PartGetDatasetInfo(file, i, name, sizeof(name), &type, &dim);
        std::cout << " Reading " << name << " " << type << " dim " << dim << std::endl;

    }

    // NOTE latest time step is number of steps - 1
    H5SetStep(file, step_num-1);


    std::vector<Scalar> X (num_particles);
    std::vector<Scalar> Y (num_particles);
    std::vector<Scalar> Z (num_particles);
    std::vector<int> id (num_particles);

    H5PartReadDataFloat32 (file, "X", &X[0]);
    H5PartReadDataFloat32 (file, "Y", &Y[0]);
    H5PartReadDataFloat32 (file, "Z", &Z[0]);

    H5PartReadDataInt32 (file, "id", &id[0]);

    auto& posf = get_objReg().create_field<VectorField>("Pos", {}, {"X", "Y", "Z"});
    auto& idf = get_objReg().create_field<IntField>("id");

    posf.reserve(posf.size());
    idf.reserve(posf.size());

    for (h5_int64_t iter=0; iter< num_particles; iter++) {
        posf.push_back({X[iter], Y[iter], Z[iter]});
        idf.push_back(id[iter]);
    }

    get_objReg().update_n_particles();

    time_graph_.set_current_timestep(step_num*write_freq_);

    H5CloseFile(file);

    // NOTE H5CloseFile seems to miss to free two pointers
    // hence the clean up here to avoid asan errors
    h5_prop_file *prop_ptr = (h5_prop_file *)prop;

    free((char *)(prop_ptr->prefix_iteration_name));
    free((h5_prop_t *)prop);
}


REGISTER_DEF_TYPE(IMPORT, HDF5Reader);
