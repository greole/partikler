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

#ifndef PARTIKLER_READERBASE_INCLUDED_H
#define PARTIKLER_READERBASE_INCLUDED_H

#include <string> // for string

#include "Models.hpp" // for TimeGraph (ptr only), Model
#include "Time.hpp"   // for TimeGraph (ptr only), Model

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

#include "H5hut.h"
#include "h5core/h5_types.h"

typedef struct h5_prop_file h5_prop_file_t;
typedef h5_prop_file_t *h5_prop_file_p;

// NOTE re-declaration of h5_prop_file to do the pointer clean-up.
//      we cant include the corresponding private/h5_types.h since it is not
//      in the include directory
struct h5_prop_file {     // file property
    h5_int64_t class_;    // property class == H5_PROP_FILE
    h5_int64_t flags;     // file access mode (read-write, readonly ...
    h5_int64_t align;     // HDF5 alignment
    h5_int64_t increment; // increment for core vfd
    h5_int64_t throttle;
#ifdef H5_HAVE_PARALLEL
    MPI_Comm comm;
#endif
    hid_t xfer_prop;             // dataset transfer properties
    hid_t access_prop;           // file access properties
    hid_t create_prop;           // file create properties
    char *prefix_iteration_name; // Prefix of step name
    int width_iteration_idx;     // pad iteration index with 0 up to this
    int flush;                   // flush iteration after writing dataset
};


class ReaderBase : public Model {

  protected:

    bool read_;

    // store reference to current TimeGraph instance
    TimeGraph &time_graph_;

  public:
    ReaderBase(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    TimeGraph &get_timeGraph() { return time_graph_; }

    bool read() {return read_;};
};

class HDF5Reader : public ReaderBase {

    REGISTER_DEC_TYPE(HDF5Reader);

private:
    std::string file_name_;

    int write_freq_;

public:
    HDF5Reader(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    void execute();
};

#endif
