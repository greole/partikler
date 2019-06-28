#ifndef FILEIO_H
#define FILEIO_H

#include <iostream>
#include <fstream>
#ifdef _WIN32
  #include<direct.h>
#else
  #include <dirent.h>
#endif
#include <string>
#include <sstream>
#include <sys/stat.h>

std::string intToStr(int number) {
    std::stringstream ss;    //create a stringstream
    ss << number;    //add number to the stream
    return ss.str();    //return a string with the contents of the stream
}

void createFolder(std::string dirname) {
    std::string str_processor= dirname;
#ifdef  _WIN32
    int stat = _mkdir((char*)str_processor.c_str());
#else  /* Annahme: Unix */
    int stat = mkdir((char*) str_processor.c_str(), 0777);
    if (stat == 0) {
        //cout<<"\n\ncreated directory "<<str_processor <<endl;
    }else{
        //cout<<"\n\nnot created directory "<<str_processor <<endl;
    }
#endif
}

// template<class T>
// void swrite(
//     const std::vector<T>& data,
//     const std::string foldername,
//     const std::string dataname,
//     const std::string type,
//     const std::string size)
// {
//     long long buf_size = data.size();
//     std::string filename = foldername + "/" + dataname + "." + type + size;
//     // std::cout<<"writing to filename: "<<filename<<std::endl;
//     std::ofstream fh;
//     fh.open(filename.c_str());
//     fh.write(reinterpret_cast<char*>(&data[0]), buf_size*sizeof(T));
//     fh.close();
// }

void prepare_data_folder (
    std::string foldername,
    int step
    ) {

    // #<{(|* pseudo time, needed for Paraview plugin|)}>#
    float realtime = (float) step;

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
}



template <class T>
void write_field (
    const std::string foldername,
    const int step,
    const std::string fieldname,
    const std::vector<T>& data,
    const std::string type,
    const std::string size
    )
{
    prepare_data_folder(foldername, step);
    std::string stepname = foldername + "/step#" + intToStr(step);

    std::vector<float> buffer(data.size());
    for (size_t i = 0; i < data.size(); i++) {
      buffer[i] = data[i];
    }

    std::cout << "writing scalar field "  << fieldname  << std::endl;
    std::string filename = stepname + "/" + fieldname + "." + type + size;

    std::ofstream fh;
    fh.open(filename.c_str());
    fh.write(reinterpret_cast<const char*>(&buffer[0]), data.size()*sizeof(T));
    fh.close();
}

template <class T>
void write_vector_field(
    const std::string foldername,
    const int step,
    const std::string fieldname,
    const std::vector<T> data,
    const std::string type,
    const std::string size,
    const std::vector<std::string> comp_names
    )
{
    prepare_data_folder(foldername, step);
    std::vector<float>buffer(data.size());
    std::string stepname = foldername + "/step#" + intToStr(step);

    size_t j = 0;
    std::cout << "writing vector field " << fieldname  << std::endl;
    for (std::string component: comp_names) {
        std::vector<float> buffer (data.size());
        for (size_t i=0; i<data.size(); i++) {
            buffer[i] = data[i][j];
        }

        std::string filename = stepname + "/" + component + ".float" + size;

        std::ofstream fh;
        fh.open(filename.c_str());
        fh.write(reinterpret_cast<const char*>(&buffer[0]), data.size()*sizeof(float));
        fh.close();
        j++;
    }
}

template <class T>
void write_point_field(
    const std::string foldername,
    const int step,
    const std::string fieldname,
    const std::vector<T> data,
    const std::string type,
    const std::string size,
    const std::vector<std::string> comp_names) {
    prepare_data_folder(foldername, step);
    std::vector<float> buffer(data.size());
    std::string stepname = foldername + "/step#" + intToStr(step);

    std::cout << "writing point field "  << fieldname  << std::endl;
    size_t j = 0;
    for (std::string component : comp_names) {
        std::vector<float> buffer(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            buffer[i] = data[i].cartesian(j);
        }

        std::string filename = stepname + "/" + component + ".float" + size;

        std::ofstream fh;
        fh.open(filename.c_str());
        fh.write(
            reinterpret_cast<const char *>(&buffer[0]),
            data.size() * sizeof(float));
        fh.close();
        j++;
    }
}

// void write_field(
//     std::string foldername,
//     int step,
//     std::string f_name,
//     std::vector<Point>& data
//     )
// {
//     prepare_data_folder(foldername, step);
//     std::vector<float>buffer(data.size());
//     std::string stepname = foldername + "/step#" + intToStr(step);
//
//     for (size_t i =0; i<data.size(); i++){
//         buffer[i]= data[i].x();
//     }
//     swrite(buffer, stepname, "X");
//
//     for (size_t i =0; i<data.size(); i++){
//         buffer[i]= data[i].y();
//     }
//     swrite(buffer, stepname, "Y");
//
//     for (size_t i =0; i<data.size(); i++){
//         buffer[i]= data[i].z();
//     }
//     swrite(buffer, stepname, "Z");
// }

#endif
