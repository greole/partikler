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

#include "Logger.hpp"
#include "Datastructures.hpp"
#include "Field.hpp"
#include "FileIO.hpp"

#include <vector>
#include <stdio.h>
#include <iostream>
#include <memory>

class ObjectRegistry {

private:

    using ObjReg = std::vector<std::unique_ptr<SPHObject>>;

    ObjReg objects_;

    size_t n_particles_ = 0;

public:

    template <class T>
    T& register_object(std::unique_ptr<SPHObject> f) {
        std::cout << "Register object " << f->get_name() << std::endl;
        int idx = objects_.size();
        objects_.push_back(std::move(f));
        return dynamic_cast<T &> ( *objects_[idx] );
    }

    // Setter

    void set_n_particles(size_t n_particles) {
        n_particles_ = n_particles;
    }

    template <class T> T &get_object(const std::string name) {
        for (auto &&f : objects_) {
            if (f->get_name() == name) {
                return dynamic_cast<T &>(*f);
            };
        }
        // for (Object *f : objects_) {
        //     if (f->get_name() == name) {
        //         return dynamic_cast<T &>(*f);
        //     };
        // }
        std::cout << "Could not find object " << name << std::endl;
    }

    ObjReg& get_objects() { return objects_; }

    bool object_exists(const std::string name) const {
        for (auto &&f : objects_) {
            if (f->get_name() == name) {
                return true;
            };
        }

        return false;
    }


    // template <class T>
    // T &get_field(const std::string name){
    //     if (object_exists(name)) return get_objects<T> (name)
    // }

    PointField &get_particle_positions() {
        return get_object<PointField &>("Pos");
    }

    // create an generic with default val
    template <class T> T &create_generic(const std::string name) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(
            std::make_unique<T>(name, GenericType, typename T::value_type()));
    }

    // create an generic with defined  val
    template <class T>
    T &create_generic(const std::string name, typename T::value_type val) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(
            std::make_unique<T>(name, GenericType, val));
    }

    template <class T>
    T &create_field(
        const std::string name,
        typename T::value_type init_value,
        const std::vector<std::string> comp_names) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_, init_value),
            comp_names, name));
    }

    template <class T>
    T &create_field(const std::string name, typename T::value_type init_value) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_, init_value),
            name));
    }

    template <class T>
    T &create_field(const std::string name) {
        if (object_exists(name)) return get_object<T>(name);
        return register_object<T>(std::make_unique<T>(
            std::vector<typename T::value_type>(n_particles_), name));
    }

    // SizeTField &create_idx_field() {
    //     // TODO move to cpp
    //     if (object_exists("idx")) return get_object<SizeTField>("idx");

    //     auto &f = create_field<SizeTField>("idx");
    //     std::cout << " n_particles_ " << n_particles_ << std::endl;

    //     for (size_t i = 0; i < n_particles_; i++) { f[i] = i; }

    //     return f;
    // }

    void write_to_disk(int timestep, std::string name="Data") {
        std::cout << "write to disk" << std::endl;
        int cur_timestep = timestep;//get_timestep();
        int write_freq = 1; // get_write_freq();
        int index_on_dist = cur_timestep / write_freq;
        std::cout << "TIMESTEP: " << cur_timestep << " " << write_freq << std::endl;

        if (cur_timestep - index_on_dist * write_freq == 0) {
            std::cout << get_objects().size() << std::endl;

            prepare_data_folder(name, index_on_dist);

            const std::string path = name + "/step#" + intToStr(index_on_dist);

            for (auto &f : get_objects()) f->write_to_disk(path);
        }
    }
};

#endif
