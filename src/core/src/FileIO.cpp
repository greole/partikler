#include "FileIO.hpp"

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

void prepare_data_folder(std::string foldername, int step) {

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
}

template <class T>
void write_field(
    const std::string foldername,
    const std::string fieldname,
    const std::vector<T> &data,
    const std::string type,
    const std::string size) {

    std::vector<float> buffer(data.size());
    for (size_t i = 0; i < data.size(); i++) buffer[i] = data[i];

    std::cout << "writing scalar field " << fieldname << data.size()<< std::endl;
    std::string filename = foldername + "/" + fieldname + "." + type + size;

    std::ofstream fh;
    fh.open(filename.c_str());
    fh.write(
        reinterpret_cast<const char *>(&buffer[0]), data.size() * sizeof(float));
    fh.close();
}

template <class T>
void write_vector_field(
    const std::string foldername,
    const std::string fieldname,
    const std::vector<T>& data,
    const std::string type,
    const std::string size,
    const std::vector<std::string> comp_names) {
    std::vector<float> buffer(data.size());

    size_t j = 0;
    std::cout << "writing vector field " << fieldname << std::endl;
    for (std::string component : comp_names) {
        std::vector<float> buffer(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            buffer[i] = data[i][j];
        }

        std::string filename = foldername + "/" + component + ".float" + size;

        std::ofstream fh;
        fh.open(filename.c_str());
        fh.write(
            reinterpret_cast<const char *>(&buffer[0]),
            data.size() * sizeof(float));
        fh.close();
        j++;
    }
}
