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

#ifndef FILEIO_H
#define FILEIO_H

#include <fstream>
#include <iostream>
#ifdef _WIN32
#include <direct.h>
#else
#include <dirent.h>
#endif
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "cgal/CGALTYPEDEFS.hpp"

std::string intToStr(int number);

void createFolder(std::string dirname);

void prepare_data_folder(std::string foldername, int step);

template <class T>
void write_scalar_field_impl(
    const T &data,
    const std::string foldername,
    const std::string fieldname,
    const std::string type,
    const std::string size) {

    std::cout << "writing scalar field " << fieldname << data.size()
              << std::endl;
    std::string filename = foldername + "/" + fieldname + "." + type + size;
    std::ofstream fh;
    fh.open(filename.c_str());
    fh.write(
        reinterpret_cast<const char *>(&data[0]),
        data.size() * sizeof(typename T::value_type));
    fh.close();
};

template <class T>
void write_vector_field_impl(
    const T &data,
    const std::string foldername,
    const std::string fieldname,
    const std::vector<std::string>& comp_names,
    const std::string type,
    const std::string size
    ) {

    size_t j = 0;
    for (std::string component : comp_names) {
        std::vector<float> buffer(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            buffer[i] = data[i][j];
        }
        write_scalar_field_impl(buffer,  foldername, component, type, size);
        j++;
    }
};

#endif
