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

#include "cgal/CGALTYPEDEFS.hpp"
#include "cgal/CGALHelper.hpp"
#include "FileIO.hpp"
#include "Logger.hpp"
#include "Helper.hpp"
#include "SearchCubes.hpp"
#include "Object.hpp"


#include <iostream>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <x86intrin.h>
#include <math.h>


static Vector zeroVec {0.,0.,0,};

// template <class T> struct NoOp {

//     NoOp(const T &val) : val_ {val} {}

//     T operator()() { return val_; }

//   private:
//     T val_;
// };

// template <class Op, class Child> struct UnaryOp {

//     UnaryOp(const Child &child) : child_(child) {}

//     auto operator()() { return op_(child_()); }

//   private:
//     Child child_;
//     Op op_;
// };

// template <class Op, class LHS, class RHS> struct BinaryOp {

//     BinaryOp(const LHS &lhs, const RHS &rhs) : lhs_(lhs), rhs_(rhs) {}

//     auto operator()() { return op_(lhs_(), rhs_()); }

//   private:
//     LHS lhs_;
//     RHS rhs_;
//     Op op_;
// };

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
    const size_t size = np.size();                                             \
    /*_Pragma("omp parallel for")*/                                            \
    for (size_t ctr = 0; ctr < size; ctr++) {                                  \
        const size_t oid = np[ctr].ownId;                                      \
        const size_t nid = np[ctr].neighId;                                    \
        /* Logger(3).check(oid <= size) << " particle oid out of bounds "; */  \
        /* Logger(3).check(nid <= size) << " particle nid out of bounds "; */  \
        op                                                                     \
    }

#define DECL_INPL_OP_SCALAR(op)                                                \
    void operator op(const ScalarField<T> &b) {                             \
        OP_FIELD_LOOP(this->f_[ctr] op b[ctr];)                                \
    }

#define DECL_INPL_OP_SCALAR_SCALAR(op)                                         \
    void operator op(const T b) { OP_LOOP(this->f_[ctr] op b;) }

#define DECL_INPL_OP_VECTOR_VECTOR(op)                          \
  void operator op(const VectorField& b) { \
    OP_FIELD_LOOP(for (int i=0;i<3;i++){this->f_[ctr][i] op b[ctr][i];} ) \
}

#define DECL_OP_SCALAR(op)                                                     \
    ScalarField<T> operator op(const ScalarField<T> &b) const {          \
        std::vector<T> ret(this->f_.size(), 0);                                \
        OP_FIELD_LOOP(ret[ctr] = this->f_[ctr] op b[ctr];)                     \
          return ScalarField<T>(ret, "tmp", this->type_);     \
    }

#define DECL_OP_VECTOR_SCALAR(op)                                              \
    VectorField operator op(const float b) const {                          \
        std::vector<Vector> ret(this->f_.size());                              \
        OP_LOOP(for (int i = 0; i < 3;                                         \
                     i++) { ret[ctr][i] = this->f_[ctr][i] op b; })            \
        return VectorField(                                                 \
            ret, {"tmpx", "tmpy", "tmpz"}, "tmp");                 \
    }

// TODO add assert here
#define DECL_OP_VECTOR_VECTOR(op)                                              \
    VectorField operator op(const VectorField &b) const {                \
        std::vector<Vector> ret(this->f_.size());                              \
        OP_LOOP(for (int i = 0; i < 3;                                         \
                     i++) { ret[ctr][i] = this->f_[ctr][i] op b[ctr][i]; })    \
        return VectorField(                                                 \
            ret, {"tmpx", "tmpy", "tmpz"}, "tmp");                 \
    }

#define DECL_OP_SCALAR_SCALAR(op)                                              \
    ScalarField<T> operator op(const T b) const {                           \
        std::vector<T> ret(this->f_.size(), 0);                                \
        OP_LOOP(ret[ctr] = this->f_[ctr] op b;)                                \
        return ScalarField<T>(ret, "tmp", this->type_);         \
    }

#define SCALAR_OPS(op)                                                         \
    DECL_INPL_OP_SCALAR(op## =)                                                \
    DECL_INPL_OP_SCALAR_SCALAR(op## =)                                         \
    DECL_OP_SCALAR_SCALAR(op)                                                  \
    DECL_OP_SCALAR(op)

// AB operations
#define OP_SCALAR_AB(name, op)                                                 \
    ScalarField<T> name##_ab(const Field<searchcubes::NeighbourPair> &np) const { \
        const size_t res_size = np.size();                              \
        std::vector<T> ret(res_size);                                          \
        OP_NEIGHBOUR_LOOP(ret[ctr] = {this->f_[oid] op this->f_[nid]};)        \
          return ScalarField<T>(ret, "tmp", this->type_);           \
    }

#define OP_VECTOR_AB(name, op)                                                 \
    VectorField name##_ab(const Field<searchcubes::NeighbourPair> &np) const { \
        const size_t ret_size = np.size();                              \
        std::vector<Vector> ret(ret_size);                                     \
        OP_NEIGHBOUR_LOOP(                                                     \
            ret[ctr] = (Vector {this->f_[oid][0] op this->f_[nid][0],          \
                                this->f_[oid][1] op this->f_[nid][1],          \
                                this->f_[oid][2] op this->f_[nid][2]});)       \
          return VectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp"); \
    }


// Move to NeighboursSearch
struct FirstLastPair{
  size_t first;
  size_t last;
};

struct VectorPair {
    // Stores the distance of particles on different
    // STL surfaces

    Vector on;
    Vector no;
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
struct Kernel {

    std::vector<float> W;

    std::vector<Vector> dWdxo;
    std::vector<Vector> dWdxn;
};


// TODO make it a member function of the fields
template <class T>
void reorder_vector(const std::vector<size_t>& idxs, std::vector<T>& vec);


//  Register a generic object
template<class T>
class Generic: public SPHObject {

protected:

    T obj_;

public:

    typedef T value_type;

    Generic(
        const std::string name,
        const std::string type,
        T obj,
        bool active=true
        ) :
        SPHObject(name, type, active),
        obj_(obj) {};

    T& operator()(){return obj_;};
};

// class SPHIOObject {
//     // convenience wrapper around SPHField classes for IO

//   public:
//     // do nothing if not implemented explicitly
//     // this allows e.g. SPHField<Kernel> to exists
//     virtual void write_to_disk(std::string path) {};
// };

template <class T> class Field : public SPHObject {

  protected:
    std::vector<T> f_;

    bool reorder_ = true;

  public:

    typedef T value_type;

    Field(
        const std::string name,
        const std::string type,
        std::vector<T> field,
        bool active = true)
        : SPHObject(name, type, active), f_(field) {};

    // a generic field
    Field(
        std::vector<T> field,
        const std::string name,
        bool active = true)
        : SPHObject(name, "generic",  active), f_(field) {};

    T &operator[](const size_t idx) { return this->f_[idx]; }
    T operator[](const size_t idx) const { return this->f_[idx]; }

    std::vector<T> &get_field() { return f_; };

    const std::vector<T> &get_field() const { return f_; };

    virtual void set_reorder(bool reorder) {

        std::cout << "set reorder_" << reorder_ << name_ << std::endl;
        reorder_ = reorder;
        std::cout << "set reorder_" << reorder_ << std::endl;
    }

    virtual void reorder(const std::vector<size_t> &idx) {
        std::cout << "reorder_" << reorder_ << std::endl;
        if (f_.size() != idx.size()) {
            logger_.warn() << " size mismatch " << name_;
        }
        if (reorder_ && f_.size() > 0) {
            logger_.info_begin() << "reordering " << name_;
            reorder_vector(idx, f_);
            logger_.info_end();
        }
    }

    void append(std::vector<T> b) {
        this->f_.insert(
            this->f_.end(),
            b.begin(),
            b.end()
            );
    };

    void set_field(std::vector<T> b) { this->f_ = b; };

    void copy_field(std::vector<T> &b) { OP_LOOP(f_[ctr] = b[ctr];) };

    void resize(size_t size) { f_.resize(size); };

    void set_uniform(T b) { OP_LOOP(f_[ctr] = b;) };

    size_t size() const { return f_.size(); };

};

// TODO make this an abstract class
template <class T> class ScalarField : public Field<T> {

  public:
    ScalarField(
        std::vector<T> field,
        const std::string name,
        const std::string type,
        bool active = true
                   )
      : Field<T>(name, type, field,  active) {};

    void set(const ScalarField &b) { OP_FIELD_LOOP(this->f_[ctr] = b[ctr];) }

    // declare all arithmetic operators
    SCALAR_OPS(+)
    SCALAR_OPS(-)
    SCALAR_OPS(*)
    SCALAR_OPS(/)

    DECL_INPL_OP_SCALAR(=)

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

    ScalarField<T> pow(const float p) const {
        OP_LOOP_RET(std::vector<float>, ret[ctr] = std::pow(this->f_[ctr], p);)
          return ScalarField<T>(ret, "tmp", this->type_);
    }

    // TODO specific class for val_ij
    void weighted_sum(
        const Field<searchcubes::NeighbourPair> &np,
        const ScalarField<float> &val_ij) {

        OP_NEIGHBOUR_LOOP(const float val = val_ij[ctr]; this->f_[oid] += val;
                          this->f_[nid] += val;)
    };
};

class IntField: public ScalarField<int> {

    public:

        IntField(
            std::vector<int> field,
            const std::string name="",
            bool active=true
                    ):
          ScalarField<int>(field, name, "int", active) {};

        void write_to_disk(std::string path);

};

class SizeTField: public ScalarField<size_t> {

    public:

        SizeTField(
                std::vector<size_t> field,
                const std::string name="",
                bool active=true
                      ):
          ScalarField<size_t>(field, name, "long",  active) {};

    void write_to_disk(std::string path);
};

class FloatField: public ScalarField<float> {

public:

  FloatField(
                std::vector<float> field,
                const std::string name = "",
                bool active=true
                ) :
    ScalarField<float>(field, name, "float",  active) {};

  FloatField(ScalarField<float> b)
    : ScalarField<float>(b.get_field(), b.get_name(), "float", true) {};

    void write_to_disk(std::string path);

};

template <class T> class ComponentField : public Field<T> {

  private:
    std::vector<std::string> comp_names_;

  public:
    ComponentField(
        std::vector<T> field,
        const std::string type_name,
        const std::string name,
        std::vector<std::string> comp_names,
        bool active=true
                      )
      : Field<T>(name, type_name, field, active),
        comp_names_(comp_names) {};

    const std::vector<std::string> get_comp_names() const { return comp_names_; };

  T& operator[](const size_t idx) { return this->f_[idx]; }

  const T operator[](const size_t idx) const { return this->f_[idx]; }


};

class VectorField : public ComponentField<Vector> {

  public:
    VectorField(
        std::vector<Vector> field,
        std::vector<std::string> comp_names,
        const std::string name = "",
        bool active=true
                   )
      : ComponentField<Vector>(field, "vector", name,  comp_names, active) {};

    void write_to_disk(std::string path);

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

  FloatField operator*(const VectorField b) const {
    // dot product
    OP_FIELD_LOOP_RET(
      std::vector<float>,
      ret[ctr] =
        this->f_[ctr][0] * b[ctr][0] +
        this->f_[ctr][1] * b[ctr][1] +
        this->f_[ctr][2] * b[ctr][2];
    )

      return FloatField(ret, "tmp");
  }

    void operator*=(const FloatField b) {
        // dot product
        OP_FIELD_LOOP(
            this->f_[ctr][0] *= b[ctr];
            this->f_[ctr][1] *= b[ctr];
            this->f_[ctr][2] *= b[ctr];
            )
    }

    void operator*=(const float b) {
        // scalar product
        OP_LOOP(
            this->f_[ctr][0] *= b;
            this->f_[ctr][1] *= b;
            this->f_[ctr][2] *= b;
            )
            }


  FloatField norm() const {
      OP_LOOP_RET(std::vector<float>,
                        ret[ctr] = std::sqrt(
                            this->f_[ctr][0] * this->f_[ctr][0] +
                            this->f_[ctr][1] * this->f_[ctr][1] +
                            this->f_[ctr][2] * this->f_[ctr][2]);)

        return FloatField(ret, "tmp");
  }


    // STL Variant of weighted_sum
  void weighted_sum(
      const Field<searchcubes::NeighbourPair> &np,
      const ScalarField<float> &val_ij,
      const Field<VectorPair> &dWdx) {

      OP_NEIGHBOUR_LOOP(for (int i = 0; i < 3; i++) {
          const float val = val_ij[ctr];
          this->f_[oid][i] += val * dWdx[ctr].on[i];
          this->f_[nid][i] += val * dWdx[ctr].no[i];
      })
  };
};

class PointField : public ComponentField<Point> {

  public:
    PointField(
        std::vector<Point> field,
        const std::string name = "",
        bool active = true
                  )
        : ComponentField<Point>(
              field,
              "Point",
              name,
              std::vector<std::string>({"X", "Y", "Z"}),
              active
                                   ) {};

    PointField(PointField &b)
        : ComponentField<Point>(
              b.get_field(),
              "Point",
              b.get_name(),
              std::vector<std::string>({"X", "Y", "Z"}),
              true) {};

    const Point &operator[](size_t id) const { return this->get_field()[id]; }

    Point &operator[](size_t id) { return this->get_field()[id]; }

    void operator+=(VectorField b) {
        const size_t size = this->f_.size();

        for (size_t i = 0; i < size; i++) {
            this->f_[i] = this->f_[i] + CGALVector(b[i][0], b[i][1], b[i][2]);
        }
  }

    void write_to_disk(std::string path);

  PointField operator-(const PointField b) const {
    // dot product
    OP_FIELD_LOOP_RET(
                      std::vector<Point>,
                      auto v = this->f_[ctr] - b[ctr];
                      ret[ctr] = Point(v[0], v[1], v[2]);
                      )

      return PointField(ret, "tmp");
  }


  FloatField norm() const {
    OP_LOOP_RET(std::vector<float>,
                ret[ctr] = std::sqrt(
                                     this->f_[ctr][0] * this->f_[ctr][0] +
                                     this->f_[ctr][1] * this->f_[ctr][1] +
                                     this->f_[ctr][2] * this->f_[ctr][2]);)

      return FloatField(ret, "tmp");
  }
};


std::ostream &operator<<(std::ostream &os, FloatField const &m);

std::ostream &operator<<(std::ostream &os, VectorField const &m);

std::ostream &operator<<(std::ostream &os, PointField const &m);

VectorField
particle_distance_vec(const PointField &pos, const Field<searchcubes::NeighbourPair> &np);

template<class T>
class EquationBase:SPHObject  {
  // Base class for equations

  // TODO Separate object registry and runTime


  private:

    T &result_;

  public:

  EquationBase(std::string name, bool active, T & result):
    SPHObject(name, "Equation", active),
    result_(result) {

    logger_.info() << " Created Equation: " << name_;
  };

    void compute();

};

#endif
