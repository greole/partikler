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

#ifndef SPHOBJECT_REGISTRY_H
#define SPHOBJECT_REGISTRY_H

#include <algorithm>          // for find_if, max
#include <ext/alloc_traits.h> // for __alloc_traits<>::value_type
#include <iostream>           // for operator<<, basic_ostream, endl, char_...
#include <map>
#include <memory> // for shared_ptr, __shared_ptr_access, make_...
#include <stdexcept>
#include <stdio.h> // for size_t
#include <string>  // for string, operator==, operator<<
#include <utility> // for move
#include <vector>  // for vector<>::iterator, vector

#include "Field.hpp"       // for PointField
#include "Logger.hpp"      // for Logger
#include "Object.hpp"      // for SPHObject, GenericType, SPHObjectType
#include "yaml-cpp/yaml.h" // for Node

class ObjectRegistry {

  private:
    using ObjReg = std::map<std::string, std::shared_ptr<SPHObject>>;

    ObjReg objects_;

    size_t n_particles_ = 0;

    Logger logger_ = {};

  public:
    template <class T> T &register_object(std::shared_ptr<SPHObject> f) {
        objects_[f->get_name()] = f;
        // return dynamic_cast<T &>(*objects_[f->get_name()]);
        // return *(objects_[f->get_name()].get());
        return dynamic_cast<T &>(*objects_[f->get_name()]);
    }

    template <class T>
    std::shared_ptr<T> register_object_get_ptr(std::shared_ptr<SPHObject> f) {
        objects_[f->get_name()] = f;
        return std::dynamic_pointer_cast<T>(objects_[f->get_name()]);
    }

    // Setter

    void set_n_particles(size_t n_particles) { n_particles_ = n_particles; }

    void update_n_particles() { n_particles_ = get_pos().size(); }

    size_t get_n_particles() {
        update_n_particles();
        return n_particles_;
    }

    template <class T> T &get_object(const std::string name) {
        // return dynamic_cast<T> (*objects_[name]);
        // return dynamic_cast<T> (objects_[name].get());
        return dynamic_cast<T &>(*objects_[name]);

        std::string error_str =
            "no object " + name + " found in object registry";
        throw std::runtime_error(error_str);
    }

    auto get_object_ptr(const std::string name) { return objects_[name]; }

    ObjReg &get_objects() { return objects_; }

    // creates a copy of the pointer of object with a new name
    void reference_clone(std::string orig_name, std::string new_name) {
        objects_[new_name] = objects_[orig_name];
    }

    bool object_exists(const std::string name) const;

    VectorField &get_pos() { return get_object<VectorField &>("Pos"); }

    // create an generic with default val
    template <class T> T &create_generic(const std::string name) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(
            std::make_unique<T>(name, GenericType, typename T::value_type()));
    }

    // create an generic with defined  val
    template <class T>
    Generic<T> &create_generic(const std::string name, T val) {
        if (object_exists(name)) return get_object<Generic<T>>(name);
        return register_object<Generic<T>>(
            std::make_unique<Generic<T>>(name, GenericType, val));
    }

    template <class T>
    T &create_field(
        const std::string name,
        typename T::value_type init_value,
        const std::vector<std::string> comp_names) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(
            std::make_unique<T>(n_particles_, init_value, name, comp_names));
    }

    template <class T>
    T &create_field(const std::string name, typename T::value_type init_value) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_, init_value),
            name));
    }

    template <class T> T &create_field(const std::string name) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_), name));
    }

    template <class T>
    T &get_or_create_model(
        const std::string model_name,
        YAML::Node parameter,
        ObjectRegistry &objReg) {
        if (object_exists(model_name)) return get_object<T>(model_name);
        return register_object<T>(
            std::make_unique<T>(model_name, parameter, objReg));
    }

    // look up velocity and create if it doesnt exist yet
    VectorField &velocity() {
        if (object_exists("u")) return get_object<VectorField>("u");
        return register_object<VectorField>(std::make_unique<VectorField>(
            n_particles_,
            zero<VectorField::value_type>::val,
            "u",
            std::vector<std::string> {"U", "V", "W"}));
    }
};

// Base class to handle materials
class Material {

  private:
    std::string name_;

    std::string type_;

    Scalar mp_;

    Scalar rho_;

  public:
    Material() {};

    Material(std::string name, std::string type, Scalar mp, Scalar rho)
        : name_(name), type_(type), mp_(mp), rho_(rho) {};

    Scalar getMp() { return mp_; };

    Scalar getRho() { return rho_; };

    std::string getName() {return name_;};
};

// Maps materials names to material objects
class MaterialMap : public SPHObject {

  private:
    std::map<std::string, Material> materials_;

  public:
    MaterialMap(const std::string name, const SPHObjectType type)
        : SPHObject(name, type) {};

    void insert(std::string material_name, Material m) {
        materials_[material_name] = m;
    };

    Material getMaterial(std::string material_name) {
        return materials_[material_name];
    };
};

// Maps field names to field ids
class FieldIdMap : public SPHObject {

  private:
    std::vector<std::string> fields_;

    std::vector<Material> material_;

  public:
    FieldIdMap(const std::string name, const SPHObjectType type)
        : SPHObject(name, type), fields_({}) {};

    // get field id from name string
    int getId(const std::string name);

    // get material from given id
    Material getMaterial(int const id) { return material_[id]; };

    int append(std::string field_name, Material m);
};

#endif
