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

#include <stdio.h> // for size_t
#include <string>  // for string
#include <vector>  // for vector

#include "Logger.hpp" // for Logger
#include "Vec3.hpp"

enum SPHObjectType {
    GenericType,
    FieldType,
    IntFieldType,
    SizeTFieldType,
    ScalarFieldType,
    VectorFieldType,
    PointFieldType,
    KernelGradientFieldType,
    EquationType,
    ModelType,
    MaterialType
};

// Template meta function to get SPHObjectType from std::vector<T>
template <class T> struct GetFieldType {
    constexpr static SPHObjectType value = GenericType;
};

template <> struct GetFieldType<std::vector<float>> {
    constexpr static SPHObjectType value = ScalarFieldType;
};

template <> struct GetFieldType<std::vector<int>> {
    constexpr static SPHObjectType value = IntFieldType;
};

template <> struct GetFieldType<std::vector<size_t>> {
    constexpr static SPHObjectType value = SizeTFieldType;
};

template <> struct GetFieldType<std::vector<Vec3>> {
    constexpr static SPHObjectType value = VectorFieldType;
};

std::string sphObjectType_to_string(SPHObjectType t);

// Base class for fields and models
class SPHObject {

  protected:
    const std::string name_;

    const SPHObjectType type_;

    //  use a copy or default here
    Logger logger_ = Logger(1);

  public:
    SPHObject() : name_(), type_(GenericType) {};

    SPHObject(const std::string name, const SPHObjectType type)
        : name_(name), type_(type) {};

    virtual ~SPHObject() = default;

    std::string get_name() const { return name_; };

    SPHObjectType get_type() const { return type_; };

    std::string get_type_str() const { return sphObjectType_to_string(type_); };
};

template <class T> class Generic : public SPHObject {

  protected:
    T obj_;

  public:
    typedef T value_type;

    Generic(const std::string name, const SPHObjectType, T obj)
        : SPHObject(name, GenericType), obj_(obj) {};

    T &operator()() { return obj_; };
};

template <class T> class EquationBase : SPHObject {
    // Base class for equations

    // TODO Separate object registry and runTime

  private:
    T &result_;

  public:
    EquationBase(std::string name, bool active, T &result)
        : SPHObject(name, EquationType), result_(result) {

        logger_.info() << " Created Equation: " << name_;
    };

    void compute();
};

#endif
