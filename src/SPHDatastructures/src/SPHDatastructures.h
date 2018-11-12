#ifndef SPHDATASTRUCTURES_H
#define SPHDATASTRUCTURES_H

#include "CGALTYPEDEFS.h"
#include "Logger.h"
#include "ParticleNeighbours.h"
#include "SearchCubes.h"

struct Kernel {

    std::vector<float> W;

    std::vector<Vector> dWdx;
};

class SPHFieldBase {


    protected:

        std::string name_;

        std::string type_;

    public:

        SPHFieldBase(std::string name, std::string type):
            name_(name),
            type_(type)
            {};

        virtual ~SPHFieldBase() = default;

        std::string get_name() const {return name_;};

        std::string get_type() const {return type_;};

};

template<class T>
class SPHField: public SPHFieldBase {

        protected:

            std::vector<T> f_;

        public:

            SPHField(
                const std::string name,
                const std::string type,
                std::vector<T> field
            ) :
                SPHFieldBase(name, type),
                f_(field) {};

            std::vector<T>& get_field() {return f_;};

            const std::vector<T>& get_field() const {return f_;};
};

class SPHIntField:public SPHField<int> {

    public:

        SPHIntField(
                 std::vector<int>& field,
                const std::string name=""):
                SPHField<int>(name, "int", field) {};
};


class SPHScalarField:public SPHField<float> {

    public:

        SPHScalarField(
            std::vector<float>& field,
            const std::string name = ""
            ) :
            SPHField<float>(name, "float", field) {};


        // Serial operations

        SPHScalarField& operator+=(SPHScalarField& b) {
            std::vector<float> & a = get_field();
            for(size_t ctr=0; ctr< a.size(); ctr++) {
                a[ctr] += b.get_field()[ctr];
            }
            return *this;
        }

        SPHScalarField operator+(SPHScalarField& b) {

            std::vector<float> ret {f_};

            for(size_t ctr=0; ctr < f_.size(); ctr++) {
                ret[ctr] += b.f_[ctr];
            }
            return SPHScalarField {ret};
        }


};

class SPHVectorField:public SPHField<Vector> {

    public:

        SPHVectorField(
                std::vector<Vector>& field,
                const std::string name=""
                ):
                SPHField<Vector>(name, "float", field) {};
};

class SPHPointField:public SPHField<Point> {

    public:

        SPHPointField(
                std::vector<Point>& field,
                const std::string name=""
                ):
                SPHField<Point>(name, "Point", field) {};
};



#endif
