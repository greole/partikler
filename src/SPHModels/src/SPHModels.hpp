#ifndef SPHMODEL_H
#define SPHMODEL_H

#include <map>
#include <stdio.h>
#include "SPHDatastructures.hpp"


#define REGISTER_DEC_TYPE(NAME)                 \
    static DerivedRegister<NAME> reg

#define REGISTER_DEF_TYPE(NAME)                 \
    DerivedRegister<NAME> NAME::reg(#NAME)


class SPHModel : public SPHObject {

  private:

  public:

    SPHModel(
        const std::string name
        ):
        SPHObject(name, "Model", false, true)
    {
        std::cout << " creating "
                  << this->get_type()
                  << " "
                  << this->get_name()
                  << std::endl;

    };

    virtual void execute() = 0;
};



struct SPHModelFactory {

    typedef std::map<std::string, SPHModel *(*)()> map_type;

  private:

    static map_type *map_;

  public:

    static SPHModel *createInstance(std::string const &model_name) {
        map_type::iterator it = getMap()->find(model_name);
        if (it == getMap()->end()) {
            std::cout << " no model named "
                      << model_name
                      << " found"
                      << std::endl;
            return 0;};
        return it->second();
    }

  protected:

    static map_type *getMap() {
        // never delete'ed. (exist until program termination)
        // because we can't guarantee correct destruction order
        if (!map_) {
            map_ = new map_type;
        }
        return map_;
    }
};

template<typename T>
SPHModel* createModel() { return new T(); }

template<typename T>
struct ModelRegister : SPHModelFactory {

    ModelRegister(std::string const& s) {
        getMap()->insert(std::make_pair(s, &createModel<T>));
    }
};

// Default Models

class SPHModelGraph : public SPHModel {

    // Defines a temporary order for submodels
     static ModelRegister<SPHModelGraph> reg;

private:

    std::vector<SPHModel*> submodels_;

public:

    SPHModelGraph(
        ) : SPHModel ("SPHModelGraph") {};

    void push_back(SPHModel &m) {

        std::cout << " Registering: "
                  << this->get_type()
                  << " "
                  << this->get_name()
                  << std::endl;

        submodels_.push_back(&m);
    };

    void execute() {
        std::cout << " Executing: "
                  << this->get_type()
                  << " "
                  << this->get_name()
                  << std::endl;

        for (auto model: submodels_) model->execute();
    };


};

#endif
