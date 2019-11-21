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
void write_int_field(
    const std::string foldername,
    const std::string fieldname,
    const std::vector<T> &data,
    const std::string type,
    const std::string size);

template <class T>
void write_field(
    const std::string foldername,
    const std::string fieldname,
    const std::vector<T> &data,
    const std::string type,
    const std::string size);

template <class T>
void write_vector_field(
    const std::string foldername,
    const std::string fieldname,
    const std::vector<T> &data,
    const std::string type,
    const std::string size,
    const std::vector<std::string> comp_names);

void write_point_field(
    const std::string & foldername,
    const std::string & fieldname,
    const std::vector<Point> & data,
    const std::string& type,
    const std::string& size,
    const std::vector<std::string>& comp_names);

#endif
