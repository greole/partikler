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

#include "HDF5Writer.hpp"

#include "Scalar.hpp"
#include "Time.hpp"

void H5PartWriteWrapper(
    h5_file_t &fh, std::string comp_name, std::vector<float> &f) {
    H5PartWriteDataFloat32(fh, comp_name.c_str(), &f[0]);
}

void H5PartWriteWrapper(
    h5_file_t &fh, std::string comp_name, std::vector<double> &f) {
    H5PartWriteDataFloat64(fh, comp_name.c_str(), &f[0]);
}

template <class T>
void HDF5Writer::write_to_disk(T const &data, h5_file_t &fh) {}

template <>
void HDF5Writer::write_to_disk(IntField const &data, h5_file_t &fh) {
    H5PartWriteDataInt32(fh, data.get_name().c_str(), &data[0]);
}

template <>
void HDF5Writer::write_to_disk(DoubleField const &data, h5_file_t &fh) {
    H5PartWriteDataFloat64(fh, data.get_name().c_str(), &data[0]);
}

template <>
void HDF5Writer::write_to_disk(FloatField const &data, h5_file_t &fh) {
    H5PartWriteDataFloat32(fh, data.get_name().c_str(), &data[0]);
}

// TODO use SFINAE here
template <>
void HDF5Writer::write_to_disk(const VectorField &data, h5_file_t &fh) {
    size_t j = 0;
    for (std::string comp : data.get_comp_names()) {
        std::vector<Scalar> buffer(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            buffer[i] = data[i][j];
        }
        H5PartWriteWrapper(fh, comp, buffer);
        j++;
    }
}

HDF5Writer::HDF5Writer(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : WriterBase(model_name, parameter, objReg),
      export_name_(parameter["name"].as<std::string>()),
      step_(0)
{}

void HDF5Writer::execute() {

    if (write()) {
        const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

        int comm_rank = 0;
        std::string fname {(export_name_ + ".h5part")};

        H5AbortOnError();
        H5SetVerbosityLevel(h5_verbosity);

        // open file and create first step
        h5_prop_t prop = H5CreateFileProp();
        H5SetPropFileCoreVFD(prop, 0);
        h5_file_t file = H5OpenFile(fname.c_str(), H5_O_RDWR, prop);
        h5_int64_t start, end;

        int cur_timestep = time_graph_.get_current_timestep();
        int index_on_dist = cur_timestep / get_write_freq();

        H5SetStep(file, index_on_dist);
        H5PartSetNumParticles(file, get_objReg().get_n_particles());

        auto &objReg = get_objReg();

        for (auto &obj : objReg.get_objects()) {

            auto name = obj.first;
            auto type = obj.second->get_type();

            std::shared_ptr<SPHObject> *obj_ptr = &obj.second;
            DISPATCH(obj_ptr, write_to_disk, type, file);
        }

        H5CloseFile(file);

        // NOTE H5CloseFile seems to miss to free two pointers
        // hence the clean up here to avoid asan errors
        h5_prop_file *prop_ptr = (h5_prop_file *)prop;

        free((char *)(prop_ptr->prefix_iteration_name));
        free((h5_prop_t *)prop);
    }
}

REGISTER_DEF_TYPE(EXPORT, HDF5Writer);
