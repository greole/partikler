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

#ifndef SPHDATASTRUCTURES_H
#define SPHDATASTRUCTURES_H

#include "CGALTYPEDEFS.hpp"
#include "CGALHelper.hpp"
#include "FileIO.hpp"
#include "Logger.hpp"
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdio.h>
#include "Helper.hpp"

#include <x86intrin.h>
#include <math.h>


static Vector zeroVec {0.,0.,0,};

// do operation op over fields f_ and b
// b must provide a size() method
// f_ and b must be active
#define OP_FIELD_LOOP(op)                                                      \
    const size_t size = this->f_.size();                                       \
    /* use logger_ instance */                                                 \
     if (this->active_ && b.active()) {                                        \
        Logger(3).check(size == b.size()) << " fields lengths not equal";      \
        _Pragma("omp parallel for")                                            \
        for (size_t ctr = 0; ctr < size; ctr++) {                              \
            op                                                                 \
        }                                                                      \
     }

// do operation op over field f_ and scalar b
// f_ must be active
#define OP_LOOP(op)                                                            \
    const size_t size = this->f_.size();                                       \
    if (this->active_) {                                                       \
    _Pragma("omp parallel for")                                                \
    for (size_t ctr = 0; ctr < size; ctr++) {                                  \
        op                                                                     \
    }                                                                          \
}

#define FOR_ALL(field, ctr)                        \
  for (size_t ctr = 0; ctr < field.size(); ctr++)  \


// like OP_LOOP but instantiates ret vector
#define OP_LOOP_RET(rettype, op)                                               \
    const size_t ret_size = this->f_.size();                                   \
    rettype ret(ret_size);                                                     \
    OP_LOOP(op)

// like OP_LOOP but instantiates ret vector
#define OP_FIELD_LOOP_RET(rettype, op)                                         \
    const size_t ret_size = this->f_.size();                                   \
    rettype ret(ret_size);                                                     \
    OP_FIELD_LOOP(op)

// do operation op over field f_ and neighbours pn
// provides ctr - counter
// oid, nid owner and neighbour ids
// size number of neighbour pairs
#define OP_NEIGHBOUR_LOOP(op)                                                  \
    const size_t size = pn.ids.size();                                       \
    /*_Pragma("omp parallel for")*/                                            \
    for (size_t ctr = 0; ctr < size; ctr++) {                                  \
        const size_t oid = pn.ids[ctr].ownId;                                      \
        const size_t nid = pn.ids[ctr].neighId;                                    \
        /* Logger(3).check(oid <= size) << " particle oid out of bounds "; */  \
        /* Logger(3).check(nid <= size) << " particle nid out of bounds "; */  \
        op                                                                     \
    }

#define DECL_INPL_OP_SCALAR(op)                                                \
    void operator op(const SPHScalarField<T> &b) {                             \
        OP_FIELD_LOOP(this->f_[ctr] op b[ctr];)                                \
    }

#define DECL_INPL_OP_SCALAR_SCALAR(op)                                         \
    void operator op(const T b) { OP_LOOP(this->f_[ctr] op b;) }

#define DECL_INPL_OP_VECTOR_VECTOR(op)                          \
  void operator op(const SPHVectorField& b) { \
    OP_FIELD_LOOP(for (int i=0;i<3;i++){this->f_[ctr][i] op b[ctr][i];} ) \
}

#define DECL_OP_SCALAR(op)                                                     \
    SPHScalarField<T> operator op(const SPHScalarField<T> &b) const {          \
        std::vector<T> ret(this->f_.size(), 0);                                \
        OP_FIELD_LOOP(ret[ctr] = this->f_[ctr] op b[ctr];)                     \
          return SPHScalarField<T>(ret, "tmp", this->type_, true);     \
    }

#define DECL_OP_VECTOR_SCALAR(op)                                              \
    SPHVectorField operator op(const float b) const {                          \
        std::vector<Vector> ret(this->f_.size());                              \
        OP_LOOP(for (int i = 0; i < 3;                                         \
                     i++) { ret[ctr][i] = this->f_[ctr][i] op b; })            \
        return SPHVectorField(                                                 \
            ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true);                 \
    }

// TODO add assert here
#define DECL_OP_VECTOR_VECTOR(op)                                              \
    SPHVectorField operator op(const SPHVectorField &b) const {                \
        std::vector<Vector> ret(this->f_.size());                              \
        OP_LOOP(for (int i = 0; i < 3;                                         \
                     i++) { ret[ctr][i] = this->f_[ctr][i] op b[ctr][i]; })    \
        return SPHVectorField(                                                 \
            ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true);                 \
    }

#define DECL_OP_SCALAR_SCALAR(op)                                              \
    SPHScalarField<T> operator op(const T b) const {                           \
        std::vector<T> ret(this->f_.size(), 0);                                \
        OP_LOOP(ret[ctr] = this->f_[ctr] op b;)                                \
        return SPHScalarField<T>(ret, "tmp", this->type_, true);         \
    }

#define SCALAR_OPS(op)                                                         \
    DECL_INPL_OP_SCALAR(op## =)                                                \
    DECL_INPL_OP_SCALAR_SCALAR(op## =)                                         \
    DECL_OP_SCALAR_SCALAR(op)                                                  \
    DECL_OP_SCALAR(op)

// AB operations
#define OP_SCALAR_AB(name, op)                                                 \
    SPHScalarField<T> name##_ab(const SortedNeighbours &pn) const {            \
        const size_t res_size = pn.ids.size();                               \
        std::vector<T> ret(res_size);                                          \
        OP_NEIGHBOUR_LOOP(ret[ctr] = {this->f_[oid] op this->f_[nid]};)        \
          return SPHScalarField<T>(ret, "tmp", this->type_, true);           \
    }

#define OP_VECTOR_AB(name, op)                                                 \
    SPHVectorField name##_ab(const SortedNeighbours &pn) const {               \
        const size_t ret_size = pn.ids.size();                               \
        std::vector<Vector> ret(ret_size);                                     \
        OP_NEIGHBOUR_LOOP(                                                     \
            ret[ctr] = (Vector {this->f_[oid][0] op this->f_[nid][0],          \
                                this->f_[oid][1] op this->f_[nid][1],          \
                                this->f_[oid][2] op this->f_[nid][2]});)       \
          return SPHVectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true); \
    }

struct Point3D {
  // 3*4bytes = 12bytes
  float x,y,z;
};

struct SubDivision {
  // 3*4bytes = 12bytes
  // Max 65535**3 cubes, which are approx (65535*3)**3 particles
  unsigned int nx, ny, nz;
};

struct FirstLastPair{
  size_t first;
  size_t last;
};

struct SearchCube {
  // A search cube stores first and last particle ids
  size_t first;
  size_t last;

  SearchCube(size_t first, size_t last) : first(first), last(last){};
};

struct SortedParticles {
  std::vector<SearchCube> searchCubes;
  std::vector<size_t> sorting_idxs;
  std::vector<Point> particles;
};

struct NeighbourPair {
    size_t ownId;
    size_t neighId;
};

struct STLSurfaceDist {
    // Stores the distance of particles on different
    // STL surfaces

    float len;
    CGALVector on;
    CGALVector no;
};

STLSurfaceDist compute_STLSurface_dist(
    Point opos, Point npos,
    Facet_handle start, Facet_handle end) {

    if (start == end) {
        // if particles on same facet return Cartesian distance
        CGALVector lenVo = npos - opos;
        CGALVector lenVn = opos - npos;
        return {(float) length(lenVo), lenVo, lenVn};
    }

    std::pair<bool, Path> ret_path = searchPath({start, end});

    if (!ret_path.first) {
        // if path search failed use Cartesian distance 
        CGALVector lenVo = npos - opos;
        CGALVector lenVn = opos - npos;
        return {(float) length(lenVo), lenVo, lenVn};
    }

    Path path = ret_path.second;

    Point nposp, oposp;
    nposp = projectedPoint(path, npos);
    reverse(path.begin(), path.end());
    oposp = projectedPoint(path, opos);

    CGALVector lenVo = nposp - opos;
    CGALVector lenVn = oposp - npos;

    return {(float) length(lenVo), lenVo, lenVn};
}


struct SortedNeighbours {
  std::vector<NeighbourPair> ids;

  // since computation of neighbouring particles
  // on different STL is expensive a vector holding
  // the STLSurfaceDist is stored here
  std::vector<STLSurfaceDist> dist;
};

struct SearchCubeDomain {

  Point3D min;
  Point3D max;

  float dx;
  float idx;

  // 32 byte

  SubDivision n;
  size_t nt;

  // 24 byte

  /* size_t padding; */

};


struct NeighbourIdHalfStencil {
  // Stores the stride of a domain, ie the difference of search cubes ids
  // for the 26 neighbour cubes. Since the stride is symmetric
  // only half of the 26 values are stored. Addition of the stencil
  // values yields the upper and subtraction the lower neighbours
  // Values are Stored as size_t for simpler addition without type conversion.

  //  Bottom        Middle        Top
  //  ^ y           ^ y           ^ y
  //  |x  x  x      |1  2  3      |10 11  12 // back
  //  |x  x  x      |x  x  0      |7   8   9 // middle
  //  |x  x  x      |x  x  x      |4   5   6 // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3
  //  Bottom        Middle        Top
  //  ^ y           ^ y           ^ y
  //  |06 07 08     | x  x  x     | x  x  x // back
  //  |03 04 05     |12  x  x     | x  x  x // middle
  //  |00 01 02     |09 10 11     | x  x  x // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3

  std::vector<size_t> stencil;

  NeighbourIdHalfStencil (size_t nx, size_t ny) {
    // TODO leave it a size_t, iterate only till 12, since
    // the stencil is symmetric anyway
    const size_t ny_nx = nx*ny;

    stencil = std::vector<size_t> {
                                   1,              // right
                                   nx - 1,         // back left
                                   nx    ,         // back centre
                                   nx + 1,         // back right
                                   ny_nx - nx - 1, // upper front left
                                   ny_nx - nx    , // upper front centre
                                   ny_nx - nx + 1, // upper front right
                                   ny_nx       -1, // upper middle left
                                   ny_nx         , // upper middle centre
                                   ny_nx      + 1, // upper middle right
                                   ny_nx + nx - 1, // upper back left
                                   ny_nx + nx    , // upper back centre
                                   ny_nx + nx + 1  // upper back right
    };
  };
};

// TODO revise
// use std::vector<Particle> as input
SearchCubeDomain
initSearchCubeDomain(const std::vector<Point> particles, float dx) {

    size_t n_particles = particles.size();
    auto bound_box = bounding_box(particles.begin(), particles.end());

    // bounds are scaled a bit to avoid particles on exact boundary 
    const float domain_extrusion = 0.01*dx; // in dx
    float min_x = bound_box.min().x() - domain_extrusion;
    float min_y = bound_box.min().y() - domain_extrusion;
    float min_z = bound_box.min().z() - domain_extrusion;
    float max_x = bound_box.max().x() + domain_extrusion;
    float max_y = bound_box.max().y() + domain_extrusion;
    float max_z = bound_box.max().z() + domain_extrusion;
    // std::cout << "initSearchCubeDomain" << " min_z " << min_z << std::endl;
    // std::cout << "initSearchCubeDomain" << " max_z " << max_z << std::endl;
    // for (auto p: particles) std::cout << p << std::endl;

    // ceil with fixed dx increases max(x,y,z) a bit
    unsigned int n_cubes_x = (unsigned int)ceil(((max_x - min_x) / dx));
    unsigned int n_cubes_y = (unsigned int)ceil(((max_y - min_y) / dx));
    unsigned int n_cubes_z = (unsigned int)ceil(((max_z - min_z) / dx));

    // Check if domain is not degenerated
    // search cubes currently assume 26 neighbour cubes
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_x >= 3)
        << " n_cubes_x less than 3";
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_y >= 3)
        << " n_cubes_y less than 3";
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_z >= 3)
        << " n_cubes_z less than 3";

    std::cout
      << " [DEBUG] SPHDatastructures.hpp::initSearchCubeDomain"
      << " min_x " << min_x
      << " min_y " << min_y
      << " min_z " << min_z
      << " max_x " << max_x
      << " max_y " << max_y
      << " max_z " << max_z
      << " n_cubes_x " << n_cubes_x
      << " n_cubes_y " << n_cubes_y
      << " n_cubes_z " << n_cubes_z
      << " dx " << dx
      << std::endl;

    const float idx = 1.0 / dx;
    return SearchCubeDomain {
        Point3D {min_x, min_y, min_z},
        Point3D {max_x, max_y, max_z},
        dx,
        idx,
        SubDivision {n_cubes_x, n_cubes_y, n_cubes_z},
        (size_t)n_cubes_x * (size_t)n_cubes_y * (size_t)n_cubes_z,
    };
}

struct Kernel {

    std::vector<float> W;

    std::vector<Vector> dWdxo;
    std::vector<Vector> dWdxn;
};


// TODO make it a member function of the fields
template <class T>
void reorder_vector(std::vector<size_t>& idxs, std::vector<T>& vec) {
  Logger logger {1};
  logger.info_begin() <<  "reorder vector";
  std::vector<T> tmp(vec.size());

  for(size_t i=0; i<idxs.size(); i++) {
    tmp[idxs[i]] = vec[i];
  }

  vec=tmp;
  logger.info_end();
}

class SPHObject {

    // Base class for fields and models

    protected:

        const std::string name_;

        const std::string type_;

        // marks an object as temporary and thus avoids
        // registration
        const bool tmp_;

        bool active_;

        //  use a copy or default here
        Logger logger_ = Logger(1);

        // TODO might need extra treatment for threads
        static std::vector<SPHObject *> objects_;

    public:

        SPHObject(
           const std::string name,
           const std::string type,
           const bool tmp=false,
           const bool active=true
                ):
            name_(name),
            type_(type),
            tmp_(tmp),
            active_(active)
            {
              // TODO dont register tmp objects
              if (!tmp_) this->register_object(*this);
            };

        virtual ~SPHObject() = default;

        // activate or deactivate object
        void set_active(bool active) {active_=active;};

        // get status
        bool active() const {return active_;};

        std::string get_name() const {return name_;};

        std::string get_type() const {return type_;};

        // register a SPHObject to object registry
        void register_object(SPHObject &f) {
          {
            // auto msg = logger_.info();
            std::cout << " Registering: "
                << f.get_type()
                << " "
                << f.get_name()
                << std::endl;
          }
          objects_.push_back(&f);
        }

        static std::vector<SPHObject *> get_objects() {
            return objects_;
        }

  template< class T>
  T& get_object(const std::string name) const {
    for (SPHObject* f: objects_) {
      if (f->get_name() == name) {return dynamic_cast<T&> (*f);};
    }
  }
};

template<class T>
class SPHField: public SPHObject {

        protected:

            std::vector<T> f_;

        public:

            SPHField(
                const std::string name,
                const std::string type,
                std::vector<T>& field,
                bool tmp=false,
                bool active=true
            ) :
              SPHObject(name, type, tmp, active),
                f_(field) {};

            std::vector<T>& get_field() {return f_;};

            const std::vector<T>& get_field() const {return f_;};

            void reorder_field(const std::vector<size_t> &idx) {
              logger_.info_begin() << "reordering " << name_;

              std::vector<T> ordered(f_.size(), f_[0]);

              for(size_t i=0; i<f_.size(); i++) {
                ordered[idx[i]] = f_[i];
              }

              f_ = ordered;
              logger_.info_end();
            }

            void set_field(std::vector<T> b) {f_ = b;};

            void set_uniform(T b) {
              OP_LOOP(f_[ctr] = b;)
            };

            size_t size() const {return f_.size();};

            void write_field(const std::string file_name) const {

                long long buf_size = f_.size();

                std::ofstream fh;
                fh.open(file_name.c_str());
                fh.write(reinterpret_cast<const char*>(&f_[0]), buf_size*sizeof(T));
                fh.close();

            }
};

// TODO make this an abstract class
template <class T> class SPHScalarField : public SPHField<T> {

  public:
    SPHScalarField(
        std::vector<T> &field,
        const std::string name,
        const std::string type,
        bool tmp = false,
        bool active = true
                   )
      : SPHField<T>(name, type, field, tmp, active) {};

    void set(const SPHScalarField &b) { OP_FIELD_LOOP(this->f_[ctr] = b[ctr];) }

    // declare all arithmetic operators
    SCALAR_OPS(+)
    SCALAR_OPS(-)
    SCALAR_OPS(*)
    SCALAR_OPS(/)

    DECL_INPL_OP_SCALAR(=)

    T &operator[](const size_t idx) { return this->f_[idx]; }
    T operator[](const size_t idx) const { return this->f_[idx]; }

    // declare all _ab operator
    OP_SCALAR_AB(add, +)
    OP_SCALAR_AB(sub, -)
    OP_SCALAR_AB(mult, *)
    OP_SCALAR_AB(div, /)

    // void set_uniform(const T b) { OP_LOOP(this->f_[ctr] = b;) }

    void lower_limit(const T b) {
        OP_LOOP(this->f_[ctr] = max(b, this->f_[ctr]);)
    }

    T get_max() {
        const size_t size = this->f_.size();
        // TODO implement data type dependent lower bound
        T ret = std::numeric_limits<T>::min();
        for (size_t ctr = 0; ctr < size; ctr++) {
            ret = max(ret, this->f_[ctr]);
        }
        return ret;
    }

    SPHScalarField<T> pow(const float p) const {
        OP_LOOP_RET(std::vector<float>, ret[ctr] = std::pow(this->f_[ctr], p);)
          return SPHScalarField<T>(ret, "tmp", this->type_, true);
    }

  // TODO specific class for val_ij
    void
    weighted_sum(const SortedNeighbours &pn, const std::vector<float> &val_ij) {

        OP_NEIGHBOUR_LOOP(const float val = val_ij[ctr]; this->f_[oid] += val;
                          this->f_[nid] += val;)
    };
};

class SPHIntField: public SPHScalarField<int> {

    public:

        SPHIntField(
            std::vector<int>& field,
            const std::string name="",
            bool tmp=false,
            bool active=true
                    ):
          SPHScalarField<int>(field, name, "int", tmp, active) {};

};

class SPHSizeTField: public SPHScalarField<size_t> {

    public:

        SPHSizeTField(
                std::vector<size_t>& field,
                const std::string name="",
                bool tmp=false,
                bool active=true
                      ):
          SPHScalarField<size_t>(field, name, "long", tmp, active) {};
};

class SPHFloatField: public SPHScalarField<float> {

public:

  SPHFloatField(
                std::vector<float>& field,
                const std::string name = "",
                bool tmp=false,
                bool active=true
                ) :
    SPHScalarField<float>(field, name, "float", tmp, active) {};

  SPHFloatField(SPHScalarField<float> b)
    : SPHScalarField<float>(b.get_field(), b.get_name(), "float", true) {};

};

template <class T> class SPHComponentField : public SPHField<T> {

  private:
    std::vector<std::string> comp_names_;

  public:
    SPHComponentField(
        std::vector<T> &field,
        const std::string type_name,
        const std::string name,
        std::vector<std::string> comp_names,
        bool tmp=false,
        bool active=true
                      )
      : SPHField<T>(name, type_name, field, tmp, active),
        comp_names_(comp_names) {};

    const std::vector<std::string> get_comp_names() const { return comp_names_; };

  T& operator[](const size_t idx) { return this->f_[idx]; }

  const T operator[](const size_t idx) const { return this->f_[idx]; }


};

class SPHVectorField : public SPHComponentField<Vector> {

  public:
    SPHVectorField(
        std::vector<Vector> &field,
        std::vector<std::string> comp_names,
        const std::string name = "",
        bool tmp=false,
        bool active=true
                   )
      : SPHComponentField<Vector>(field, "vector", name,  comp_names, tmp, active) {};

    // TODO remove it and replace by set for ABC
    void set(Vector b) { OP_LOOP(f_[ctr] = b;) };

  //   void set(const SPHVectorField &b) {
  //     OP_FIELD_LOOP(
  //           this->f_[ctr][0] = b[ctr][0];
  //           this->f_[ctr][1] = b[ctr][1];
  //           this->f_[ctr][2] = b[ctr][2];
  //           )
  // }

  // vector scalar operations
  DECL_OP_VECTOR_SCALAR(*)
  DECL_OP_VECTOR_SCALAR(/)

  // vector vector operations
  DECL_OP_VECTOR_VECTOR(+)
  DECL_OP_VECTOR_VECTOR(-)

  DECL_INPL_OP_VECTOR_VECTOR(+=)
  DECL_INPL_OP_VECTOR_VECTOR(-=)
  DECL_INPL_OP_VECTOR_VECTOR(=)

// declare assignment operator
  // void operator=(const SPHVectorField &b) {
  // std::vector<Vector> ret(this->f_.size());
  // OP_LOOP(for (int i = 0; i < 3; i++) { this->f_[ctr][i] = b[ctr][i]; })
  //   }


  OP_VECTOR_AB(sub, -)
  OP_VECTOR_AB(add, +)
  OP_VECTOR_AB(mult, *)
  OP_VECTOR_AB(div, /)

  SPHFloatField operator*(const SPHVectorField b) const {
    // dot product
    OP_FIELD_LOOP_RET(
      std::vector<float>,
      ret[ctr] =
        this->f_[ctr][0] * b[ctr][0] +
        this->f_[ctr][1] * b[ctr][1] +
        this->f_[ctr][2] * b[ctr][2];
    )

      return SPHFloatField(ret, "tmp", true);
  }

    void operator*=(const SPHFloatField b) {
        // dot product
        OP_FIELD_LOOP(
            this->f_[ctr][0] *= b[ctr];
            this->f_[ctr][1] *= b[ctr];
            this->f_[ctr][2] *= b[ctr];
            )
    }

  SPHFloatField norm() const {
      OP_LOOP_RET(std::vector<float>,
                        ret[ctr] = std::sqrt(
                            this->f_[ctr][0] * this->f_[ctr][0] +
                            this->f_[ctr][1] * this->f_[ctr][1] +
                            this->f_[ctr][2] * this->f_[ctr][2]);)

        return SPHFloatField(ret, "tmp", true);
  }
};

class SPHPointField : public SPHComponentField<Point> {

  public:
    SPHPointField(
        std::vector<Point> &field,
        const std::string name = "",
        bool tmp = false,
        bool active = true
                  )
        : SPHComponentField<Point>(
              field,
              "Point",
              name,
              std::vector<std::string>({"X", "Y", "Z"}),
              tmp,
              active
                                   ) {};

    SPHPointField(SPHPointField &b)
        : SPHComponentField<Point>(
              b.get_field(),
              "Point",
              b.get_name(),
              std::vector<std::string>({"X", "Y", "Z"}),
              true) {};

    const Point &operator[](size_t id) const { return this->get_field()[id]; }

    Point &operator[](size_t id) { return this->get_field()[id]; }

    void operator+=(SPHVectorField b) {
        const size_t size = this->f_.size();

        for (size_t i = 0; i < size; i++) {
            this->f_[i] = this->f_[i] + CGALVector(b[i][0], b[i][1], b[i][2]);
        }
  }

  SPHPointField operator-(const SPHPointField b) const {
    // dot product
    OP_FIELD_LOOP_RET(
                      std::vector<Point>,
                      auto v = this->f_[ctr] - b[ctr];
                      ret[ctr] = Point(v[0], v[1], v[2]);
                      )

      return SPHPointField(ret, "tmp", true);
  }


  SPHFloatField norm() const {
    OP_LOOP_RET(std::vector<float>,
                ret[ctr] = std::sqrt(
                                     this->f_[ctr][0] * this->f_[ctr][0] +
                                     this->f_[ctr][1] * this->f_[ctr][1] +
                                     this->f_[ctr][2] * this->f_[ctr][2]);)

      return SPHFloatField(ret, "tmp", true);
  }
};


std::ostream &operator<<(std::ostream &os, SPHFloatField const &m) {

    auto fs = m.get_field();
    os
     << "\n"
     << m.get_name() << " " << fs.size()
     << std::endl;

    size_t id = 0;
    for (auto f : fs) {
      os << m.get_name() << "Particle " << id << " " << f << "\n" << std::flush;
      id++;
    }
    return os;

}

std::ostream &operator<<(std::ostream &os, SPHVectorField const &m) {

  auto fs = m.get_field();
  os 
     << "\n"
     << m.get_name() << " " << fs.size()
     << std::endl;
  size_t id = 0;
  for (auto f : fs) {os
      << m.get_name()
      << " Particle " << id << " ("
      << f[0] << " "
      << f[1] << " "
      << f[2] << ")\n" << std::flush;
    id++;
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, SPHPointField const &m) {

  auto fs = m.get_field();
  os 
     << "\n"
     << m.get_name() << " " << fs.size()
     << std::endl;
  size_t id = 0;


  for (auto f : fs) {os
                      << m.get_name()
                      << " Particle " << id << " ("
                      << f[0] << " "
                      << f[1] << " "
                      << f[2] << ")\n" << std::flush;
      id++;
    }
  return os;
}

SPHVectorField
particle_distance_vec(const SPHPointField &pos, const SortedNeighbours &pn) {

    const size_t ret_size {pn.ids.size()};

    std::vector<Vector> ret(ret_size);

    OP_NEIGHBOUR_LOOP(const Point &opos = pos[oid];
                      const Point &npos = pos[nid];
                      const CGALVector lenV = npos - opos;

                      ret[ctr][0] = lenV[0];
                      ret[ctr][1] = lenV[1];
                      ret[ctr][2] = lenV[2];)

      return SPHVectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true);
}



template<class T>
class SPHEquationBase:SPHObject  {
  // Base class for equations

  // TODO Separate object registry and runTime


  private:

    T &result_;

  public:

  SPHEquationBase(std::string name, bool active, T & result):
    SPHObject(name, "SPHEquation", active),
    result_(result) {

    logger_.info() << " Created Equation: " << name_;
  };

    void compute();

};
#endif
