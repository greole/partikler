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

#ifndef PARTIKLER_WRITER_INCLUDED_H
#define PARTIKLER_WRITER_INCLUDED_H

#include <string> // for string

#include "Models.hpp" // for TimeGraph (ptr only), Model
#include "Time.hpp"   // for TimeGraph (ptr only), Model

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class WriterBase : public Model {

  protected:
    int write_freq_;

    int last_write_;

    // store reference to current TimeGraph instance
    TimeGraph &time_graph_;

  public:
    WriterBase(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    TimeGraph &get_timeGraph() { return time_graph_; }

    int get_write_freq() { return write_freq_; };

    bool write();
};

#endif
