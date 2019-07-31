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

class SPHObjectRegistry {

private:

    using ObjReg = std::vector<std::unique_ptr<SPHObject>>;

    ObjReg objects_;

public:

    template <class T>
    T& register_object(std::unique_ptr<SPHObject> f) {
        std::cout << "Register object " << f->get_name() << std::endl;
        int idx = objects_.size();
        objects_.push_back(std::move(f));
        return dynamic_cast<T &> ( *objects_[idx] );
    }


    // void register_object(std::unique_ptr<SPHObject> f) {
    //     objects_.push_back(std::move(f));
    // }

    template <class T> T &get_object(const std::string name) const {
        for (auto &&f : objects_) {
            if (f->get_name() == name) {
                return dynamic_cast<T &>(*f);
            };
        }
        // for (SPHObject *f : objects_) {
        //     if (f->get_name() == name) {
        //         return dynamic_cast<T &>(*f);
        //     };
        // }
        std::cout << "Could not find object " << name << std::endl;
    }

    ObjReg& get_objects() { return objects_; }

};

#endif
