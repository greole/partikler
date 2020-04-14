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

#include "SuperSPHWriter.hpp"

#include "Scalar.hpp"
#include "Time.hpp"

std::string intToStr(int number) {
    std::stringstream ss; // create a stringstream
    ss << number;         // add number to the stream
    return ss.str();      // return a string with the contents of the stream
}

void createFolder(std::string dirname) {
    std::string str_processor = dirname;
#ifdef _WIN32
    int stat = _mkdir((char *)str_processor.c_str());
#else /* Annahme: Unix */
    int stat = mkdir((char *)str_processor.c_str(), 0777);
    if (stat == 0) {
        // cout<<"\n\ncreated directory "<<str_processor <<endl;
    } else {
        // cout<<"\n\nnot created directory "<<str_processor <<endl;
    }
#endif
}

std::string prepare_data_folder(std::string foldername, int step) {

    // #<{(|* pseudo time, needed for Paraview plugin|)}>#
    float realtime = (float)step;

    // #<{(|* main result folder |)}>#
    // //std::cout<<"foldername: "<< foldername<<std::endl;
    createFolder(foldername);

    // #<{(|* step folders: #0, #1, #2, ... |)}>#
    std::string stepname = foldername + "/step#" + intToStr(step);
    createFolder(stepname);

    // #<{(|* times.txt: Time in ascii format|)}>#
    std::ofstream ftimes;
    std::string nametimes = stepname + "/times.txt";
    ftimes.open(nametimes.c_str());
    ftimes << realtime;
    ftimes.close();

    // #<{(|* .sph file: Format specifier for Paraview |)}>#
    std::ofstream fdotfile;
    std::string dotName = foldername + ".sph";
    fdotfile.open(dotName.c_str(), std::fstream::out);
    fdotfile.close();

    return stepname;
}

std::string field_type_to_str(SPHObjectType t) {
    switch (t) {
    case (IntFieldType):
        return "int";
    case (SizeTFieldType):
        return "int";
    // case (ScalarFieldType):
    //     return "float";
    // case (PointFieldType):
    //     return "float";
    // case (VectorFieldType):
    //     return "float";
    default:
        return "float";
    }
}

template <class T>
void write_to_disk_impl(
    T const &data,
    const std::string path,
    const std::string name,
    const SPHObjectType type) {

    std::string filename =
        path + "/" + name + "." + field_type_to_str(type) + "32";
    std::ofstream fh;
    fh.open(filename.c_str());
    fh.write(
        reinterpret_cast<const char *>(&data[0]),
        data.size() * sizeof(typename T::value_type));
    fh.close();
}

template <class T>
void SuperSPHWriter::write_to_disk(T const &data, const std::string path) {
    write_to_disk_impl(data, path, data.get_name(), data.get_type());
}

// Dont write size_t fields since the SuperSPHWriter cant read them
// template <>
// void SuperSPHWriter::write_to_disk(T const &data, const std::string path) {
//     write_to_disk_impl(data, path, data.get_name(), data.get_type());
// }

// template <>
// void SuperSPHWriter::write_to_disk<PointField>(
//     const PointField &data, const std::string path) {}

template <>
void SuperSPHWriter::write_to_disk<SizeTField>(
    const SizeTField &data, const std::string path) {

    std::vector<int> buffer(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        buffer[i] = (int)data[i];
    }
    write_to_disk_impl(buffer, path, data.get_name(), data.get_type());
}

template <>
void SuperSPHWriter::write_to_disk<VectorField>(
    const VectorField &data, const std::string path) {

    size_t j = 0;
    for (std::string comp : data.get_comp_names()) {
        std::vector<Scalar> buffer(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            buffer[i] = data[i][j];
        }
        write_to_disk_impl(buffer, path, comp, data.get_type());
        j++;
    }
}

SuperSPHWriter::SuperSPHWriter(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : WriterBase(model_name, parameter, objReg),
      export_name_(parameter["name"].as<std::string>()) {}

void SuperSPHWriter::execute() {

    if (write()) {
        int cur_timestep = get_timeGraph().get_current_timestep();
        int index_on_dist = cur_timestep / get_write_freq();

        auto stepname = prepare_data_folder(export_name_, index_on_dist);

        auto &objReg = get_objReg();

        for (auto &obj : objReg.get_objects()) {

            auto name = obj.first;
            auto type = obj.second->get_type();

            // TODO
            std::shared_ptr<SPHObject> *obj_ptr = &obj.second;
            DISPATCH(obj_ptr, write_to_disk, type, stepname);
        }
    }
}

REGISTER_DEF_TYPE(EXPORT, SuperSPHWriter);
