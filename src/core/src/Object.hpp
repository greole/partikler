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

#ifndef SPHOBJECT_H
#define SPHOBJECT_H

#include "Logger.hpp"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <memory>

class SPHObject {

    // Base class for fields and models

  protected:

    const std::string name_;

    const std::string type_;

    // marks an object as temporary and thus avoids
    // registration
    // const bool tmp_;

    bool active_;

    //  use a copy or default here
    Logger logger_ = Logger(1);

    static std::vector<SPHObject*> objects_;

  public:
    SPHObject(
        const std::string name,
        const std::string type,
        const bool active = true)
        : name_(name), type_(type), active_(active) {
        // TODO dont register tmp objects
        // if (!tmp_) this->register_object(*this);
    };

    virtual ~SPHObject() = default;

    // activate or deactivate object
    void set_active(bool active) { active_ = active; };

    // get status
    bool active() const { return active_; };

    std::string get_name() const { return name_; };

    std::string get_type() const { return type_; };


    // reorder after particle sorting
    virtual void reorder(const std::vector<size_t> &idxs) {};

    virtual void write_to_disk(std::string path) {};

};

#endif
