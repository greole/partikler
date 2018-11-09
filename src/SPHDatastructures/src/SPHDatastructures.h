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


    private:

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

        private:

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
};

class SPHIntField:public SPHField<int> {

    public:

        SPHIntField(
                 std::vector<int>& field,
                const std::string name=""):
                SPHField<int>(name, "Int", field) {};
};


class SPHScalarField:public SPHField<float> {

    public:

        SPHScalarField(
            std::vector<float>& field,
            const std::string name = ""
            ) :
            SPHField<float>(name, "Scalar", field) {};


        // Serial operations

        SPHScalarField& operator+=(SPHScalarField& b) {
            std::vector<float> & a = get_field();
            for(size_t ctr=0; ctr< a.size(); ctr++) {
                a[ctr] += b.get_field()[ctr];
            }
            return *this;
        }

        SPHScalarField operator+(SPHScalarField& b) {
            std::vector<float>& af =   get_field();
            std::vector<float>& bf = b.get_field();
            std::vector<float> ret (af);

            for(size_t ctr=0; ctr < af.size(); ctr++) {
                ret[ctr] += bf[ctr];
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
                SPHField<Vector>(name, "Vector", field) {};
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
