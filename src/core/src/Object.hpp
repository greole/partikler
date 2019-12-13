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

// Base class for fields and models
class SPHObject {


  protected:

    const std::string name_;

    const std::string type_;

    //  use a copy or default here
    Logger logger_ = Logger(1);

  public:

    SPHObject() {};

    SPHObject(
        const std::string name,
        const std::string type
        )
        : name_(name), type_(type){
    };

    virtual ~SPHObject() = default;

    std::string get_name() const { return name_; };

    std::string get_type() const { return type_; };

    // reorder after particle sorting
    virtual void reorder(const std::vector<size_t> &idxs) {};

    virtual void write_to_disk(std::string path) {};

};

template<class T>
class Generic: public SPHObject {

protected:

    T obj_;

public:

    typedef T value_type;

    Generic(
        const std::string name,
        const std::string type,
        T obj
        ) :
        SPHObject(name, type),
        obj_(obj) {};

    T& operator()(){return obj_;};
};


template<class T>
class EquationBase:SPHObject  {
    // Base class for equations

    // TODO Separate object registry and runTime


private:

    T &result_;

public:

    EquationBase(std::string name, bool active, T & result):
        SPHObject(name, "Equation"),
        result_(result) {

        logger_.info() << " Created Equation: " << name_;
    };

    void compute();

};

#endif
