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

#include "CGALTYPEDEFS.hpp"

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
