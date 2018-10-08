#include <CGAL/Simple_cartesian.h>

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

typedef CGAL::Simple_cartesian<double>    K;
typedef K::Point_3                        Point;

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

void swrite(
    std::vector<float>& data,
    std::string foldername,
    std::string dataname)
{
    long long buf_size = data.size();
    std::string filename = foldername + "/" + dataname + ".float32";
    std::cout<<"writing to filename: "<<filename<<std::endl;
    std::ofstream fh;
    fh.open(filename.c_str());
    fh.write(reinterpret_cast<char*>(&data[0]), buf_size*sizeof(float));
    fh.close();
}

void writeData_SPH(
    std::string foldername,
    int step,
    std::vector<Point> data)
{
    /** pseudo time, needed for Paraview plugin*/
    float realtime = (float) step;

    /** main result folder */
    //std::cout<<"foldername: "<< foldername<<std::endl;
    createFolder(foldername);

    /** step folders: #0, #1, #2, ... */
    std::string stepname = foldername + "/step#" + intToStr(step);
    createFolder(stepname);

    /** times.txt: Time in ascii format*/
    std::ofstream ftimes;
    std::string nametimes = stepname + "/times.txt";
    ftimes.open(nametimes.c_str());
    ftimes << realtime;
    ftimes.close();

    /** .sph file: Format specifier for Paraview */
    std::ofstream fdotfile;
    std::string dotName = foldername + ".sph";
    fdotfile.open(dotName.c_str(), std::fstream::out);
    fdotfile.close();

    std::cout<<"data.size(): "<<data.size()<<std::endl;

    /** fill buffer + write */
    std::vector<float>buffer(data.size());
    for (size_t i =0; i<data.size(); i++){
        buffer[i]= data[i].x();
    }
    swrite(buffer, stepname, "X");

    for (size_t i =0; i<data.size(); i++){
        buffer[i]= data[i].y();
    }
    swrite(buffer, stepname, "Y");

    for (size_t i =0; i<data.size(); i++){
        buffer[i]= data[i].z();
    }
    swrite(buffer, stepname, "Z");
}

