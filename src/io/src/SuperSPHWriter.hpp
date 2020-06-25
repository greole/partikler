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

#ifndef PARTIKLER_SUPERSPHWRITER_INCLUDED_H
#define PARTIKLER_SUPERSPHWRITER_INCLUDED_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string> // for string
#include <sys/stat.h>

#include "Field.hpp"
#include "Models.hpp"     // for ModelRegister (ptr only), REGISTER_DEC_TYPE
#include "WriterBase.hpp" // for WriterBase
#include "yaml-cpp/yaml.h"

class ObjectRegistry;
namespace YAML {
class Node;
} // namespace YAML

class SuperSPHWriter : public WriterBase {

    REGISTER_DEC_TYPE(SuperSPHWriter);

  private:
    std::string export_name_;

  public:
    SuperSPHWriter(
        const std::string &model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg);

    template <class T> void write_to_disk(T const &t, const std::string path);

    void execute();
};

#endif
