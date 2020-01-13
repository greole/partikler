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

#ifndef PARTIKLER_HDF5_INCLUDED_H
#define PARTIKLER_HDF5_INCLUDED_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string> // for string
#include <sys/stat.h>

#include "H5hut.h"
#include "h5core/h5_types.h"

#include "Field.hpp"
#include "Models.hpp"     // for ModelRegister (ptr only), REGISTER_DEC_TYPE
#include "WriterBase.hpp" // for WriterBase
#include "cgal/CGALHelper.hpp"
#include "yaml-cpp/yaml.h"

// NOTE re-declaration of h5_prop_file to do the pointer clean-up.
//      we cant include the corresponding private/h5_types.h since it is not
//      in the include directory
struct h5_prop_file {                  // file property
        h5_int64_t class_;             // property class == H5_PROP_FILE
        h5_int64_t flags;              // file access mode (read-write, readonly ...
        h5_int64_t align;              // HDF5 alignment
       h5_int64_t increment;           // increment for core vfd
        h5_int64_t throttle;
#ifdef H5_HAVE_PARALLEL
        MPI_Comm comm;
#endif
       hid_t   xfer_prop;              // dataset transfer properties
       hid_t   access_prop;            // file access properties
       hid_t   create_prop;            // file create properties
       char*   prefix_iteration_name;  // Prefix of step name
       int     width_iteration_idx;    // pad iteration index with 0 up to this
       int     flush;                  // flush iteration after writing dataset
};

typedef struct h5_prop_file h5_prop_file_t;
typedef h5_prop_file_t* h5_prop_file_p;


class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class HDF5Writer : public WriterBase {

    REGISTER_DEC_TYPE(HDF5Writer);

  private:
    std::string export_name_;

  public:
    HDF5Writer(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    template <class T> void write_to_disk(T const &t, h5_file_t& fh);

    void execute();
};

#endif
