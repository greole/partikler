#ifndef SPHDATASTRUCTURES_H
#define SPHDATASTRUCTURES_H

#include "CGALTYPEDEFS.hpp"
#include "FileIO.hpp"
#include "Logger.hpp"
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdio.h>
// #include "ParticleNeighbours.hpp"
// #include "SearchCubes.h"
// #include "SearchCubes2.h"
// #include "SearchCubes3.h"
// #include "SearchCubes4.h"
#include "Helper.hpp"

#include <x86intrin.h>
#include <math.h>

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

struct SortedNeighbours {
  // Each particle id has a first and last neighbour
  // stored in firstLast, here SearchCube5 is reused
  /* std::vector<packed> ownId; */
  std::vector<size_t> ownId;
  std::vector<size_t> neighId;
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
    float min_x = bound_box.min().x() - 0.01*dx;
    float min_y = bound_box.min().y() - 0.01*dx;
    float min_z = bound_box.min().z() - 0.01*dx;
    float max_x = bound_box.max().x() + 0.01*dx;
    float max_y = bound_box.max().y() + 0.01*dx;
    float max_z = bound_box.max().z() + 0.01*dx;
    // std::cout << "initSearchCubeDomain" << " min_z " << min_z << std::endl;
    // std::cout << "initSearchCubeDomain" << " max_z " << max_z << std::endl;
    // for (auto p: particles) std::cout << p << std::endl;

    // ceil with fixed dx increases max(x,y,z) a bit
    unsigned int n_cubes_x = (unsigned int)ceil(((max_x - min_x) / dx));
    unsigned int n_cubes_y = (unsigned int)ceil(((max_y - min_y) / dx));
    unsigned int n_cubes_z = (unsigned int)ceil(((max_z - min_z) / dx));

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

    std::vector<Vector> dWdx;
};


template <class T>
void reorder_vector(std::vector<size_t>& idxs, std::vector<T>& vec) {
  std::cout << "[DEBUG] SPHDataStructures.hpp reorder_vector" << std::endl;
  std::vector<T> tmp(vec.size());

  for(size_t i=0; i<idxs.size(); i++) {
    tmp[idxs[i]] = vec[i];
  }

  vec=tmp;
  std::cout << "[DEBUG] SPHDataStructures.hpp reorder_vector done" << std::endl;
}


class SPHFieldBase {


    protected:

        const std::string name_;

        const std::string type_;

    public:

        SPHFieldBase(
           const std::string name,
           const std::string type
                ):
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
                std::vector<T>& field
            ) :
                SPHFieldBase(name, type),
                f_(field) {};

            std::vector<T>& get_field() {return f_;};

            const std::vector<T>& get_field() const {return f_;};

            void reorder_field(const std::vector<size_t> idx) {

              std::vector<T> ordered(f_.size(), f_[0]);

              for(size_t i=0; i<f_.size(); i++) {
                ordered[idx[i]] = f_[i];
              }

              f_ = ordered;
            }

            void set_field(std::vector<T> b) {f_ = b;};

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
template <class T>
class SPHScalarField: public SPHField<T> {


    public:

        SPHScalarField(
            std::vector<T>& field,
            const std::string name,
            const std::string type
            ):
            SPHField<T>(name, type, field) {};

        void set(const SPHScalarField b) {
            const std::vector<T> &bf = b.get_field();
            const size_t size = this->f_.size();
            for (size_t ctr = 0; ctr < size; ctr++) {
                this->f_[ctr] = bf[ctr];
            }
        }

        void operator+=(const SPHScalarField<T> b) {
            const std::vector<T> & bf = b.get_field();
            const size_t size = this->f_.size();
            for(size_t ctr=0; ctr < size;  ctr++) {
                this->f_[ctr] += bf[ctr];
            }
        }

        void operator*=(const SPHScalarField<T> b) {
            const std::vector<T> & bf = b.get_field();
            const size_t size = this->f_.size();

            for(size_t ctr=0; ctr < size; ctr++) {
                this->f_[ctr+0] *= bf[ctr+0];
            }

        }

        void operator*=(const T b) {
            const size_t size = this->f_.size();
            for(size_t ctr=0; ctr < size; ctr++) {
                this->f_[ctr] *= b;
            }

        }

        void set_uniform(const T b) {
          const size_t size = this->f_.size();
          for(size_t ctr=0; ctr < size; ctr++) {
            this->f_[ctr] = b;
          }

        }

      void lower_limit(const T b) {
        const size_t size = this->f_.size();
        for(size_t ctr=0; ctr < size; ctr++) {
          this->f_[ctr] = max(b, this->f_[ctr] );
        }

      }

        // T operator[](const size_t idx) { return this->f_[idx]; }

        T& operator[](const size_t idx) { return this->f_[idx]; }

        SPHScalarField<T> operator/(const SPHScalarField<T> b) const {
            std::vector<T> ret(this->f_.size(), 0);
            const std::vector<T> &bf = b.get_field();

            const size_t size = this->f_.size();

            for (size_t i = 0; i < size; i++) {
                ret[i] = this->f_[i] / bf[i];
            }

            return SPHScalarField<T>(ret, "tmp", this->type_);
        }

        SPHScalarField<T> operator/(const T b) const {
            std::vector<T> ret(this->f_.size(), 0);

            const size_t size = this->f_.size();

            for (size_t i = 0; i < size; i++) {
                ret[i] = this->f_[i] / b;
            }

            return SPHScalarField<T>(ret, "tmp", this->type_);
        }

        SPHScalarField<T> operator-(const T b) const {
            std::vector<T> ret(this->f_.size(), 0);

            const size_t size = this->f_.size();

            for (size_t i = 0; i < size; i++) {
                ret[i] = this->f_[i] - b;
            }

            return SPHScalarField<T>(ret, "tmp", this->type_);
        }

        SPHScalarField<T> operator+(const T b) const {
            std::vector<T> ret(this->f_.size(), 0);

            const size_t size = this->f_.size();

            for (size_t i = 0; i < size; i++) {
                ret[i] = this->f_[i] + b;
            }

            return SPHScalarField<T>(ret, "tmp", this->type_);
        }

        SPHScalarField<T> operator*(const T b) const {
            std::vector<T> ret(this->f_.size(), 0);

            const size_t size = this->f_.size();

            for (size_t i = 0; i < size; i++) {
                ret[i] = this->f_[i] * b;
            }

            return SPHScalarField<T>(ret, "tmp", this->type_);
        }

  SPHScalarField<T> operator*(const SPHScalarField<T> b) const {
    std::vector<T> ret(this->f_.size(), 0);
    const std::vector<T> &bf = b.get_field();

    const size_t size = this->f_.size();

    for (size_t i = 0; i < size; i++) {
      ret[i] = this->f_[i] * bf[i];
    }

    return SPHScalarField<T>(ret, "tmp", this->type_);
  }


        SPHScalarField<T> pow(const float p) const {
            std::vector<T> ret(this->f_.size(), 0);

            const size_t size = this->f_.size();

            for (size_t i = 0; i < size; i++) {
                ret[i] = std::pow(this->f_[i], p);
            }

            return SPHScalarField<T>(ret, "tmp", this->type_);
        }

//         void simd_inplace_mult(const SPHScalarField<T> &b) {
//             const std::vector<T> &bf = b.get_field();
//             const size_t size = this->f_.size();
//             for (size_t ctr = 0; ctr < size / 8 * 8; ctr += 8) {
//             // __m256 as = _mm256_load_ps[&(this->f_[ctr])];
//             // __m256 bs = _mm256_load_ps[&(bf[ctr])];
//             #ifdef SIMD

//             __m256 as = _mm256_setr_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0);
//             __m256 bs = _mm256_setr_ps(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0);
//             __m256 result = _mm256_mul_ps(as, bs);

//             float* res = (float*)&result;
//             this->f_[ctr+0] = res[0];
//             this->f_[ctr+1] = res[1];
//             this->f_[ctr+2] = res[2];
//             this->f_[ctr+3] = res[3];
//             this->f_[ctr+4] = res[4];
//             this->f_[ctr+5] = res[5];
//             this->f_[ctr+6] = res[6];
//             this->f_[ctr+7] = res[7];
// #else
// #endif
//         }
//     };

    // SPHScalarField<T>
    // op_ab(const SortedNeighbours &particle_neighbours, SPHScalarField<T> b,
    // functor) const {};

    SPHScalarField<T>
    add_ab(const SortedNeighbours &particle_neighbours) const {

        const size_t res_size = particle_neighbours.ownId.size();

        std::vector<T> ret(res_size);

        for (size_t ctr = 0; ctr < res_size; ctr++) {
            const size_t oid = particle_neighbours.ownId[ctr];
            const size_t nid = particle_neighbours.neighId[ctr];
            ret[ctr] = {this->f_[oid] + this->f_[nid]};
        }

        return SPHScalarField<T>(ret, "tmp", this->type_);
    }

  SPHScalarField<T>
  sub_ab(const SortedNeighbours particle_neighbours) const {

    const size_t res_size = particle_neighbours.ownId.size();

    std::vector<T> ret(res_size);

    for (size_t ctr = 0; ctr < res_size; ctr++) {
      const size_t oid = particle_neighbours.ownId[ctr];
      const size_t nid = particle_neighbours.neighId[ctr];
      ret[ctr] = {this->f_[oid] + this->f_[nid]};
    }

    return SPHScalarField<T>(ret, "tmp", this->type_);
  }

  SPHScalarField<T>
  mult_ab(const SortedNeighbours &particle_neighbours) const {

    const size_t res_size = particle_neighbours.ownId.size();

    std::vector<T> ret(res_size);

    for (size_t ctr = 0; ctr < res_size; ctr++) {
      const size_t oid = particle_neighbours.ownId[ctr];
      const size_t nid = particle_neighbours.neighId[ctr];
      ret[ctr] = {this->f_[oid] * this->f_[nid]};
    }

    return SPHScalarField<T>(ret, "tmp", this->type_);
  }

  // TODO specific class for val_ij
  void weighted_sum(
      const SortedNeighbours &neighbours, const std::vector<float> &val_ij) {

      for (size_t i = 0; i < neighbours.ownId.size(); i++) {
          size_t ownId = neighbours.ownId[i];
          size_t neighId = neighbours.neighId[i];

          const float val = val_ij[i];
          // if (val == 0.0) {std::cout << "zero kernel " << ownId << " " <<
          // neighId << std::endl;}
          this->f_[ownId]   += val;
          this->f_[neighId] += val;
      }

      //
      size_t ctr = 0;
      for (auto f : this->f_) {

          if (f == 0) {
              std::cout
                << "SPHDataStructures:487  Particle Id: "
                << ctr <<  " weigthed sum==0" << std::endl;
          }
          ctr++;
      }
  };

    // void weighted_sum(
    //     const SortedNeighbours &neighbours,
    //     SPHScalarField<T> &val_ij,
    //     Kernel &kernel);

    // void weighted_sum(
    //     const SortedNeighbours &neighbours,
    //     SPHScalarField<T> &val_ij,
    //     Kernel &kernel) {

    //     for (size_t i = 0; i < neighbours.ownId.size(); i++) {
    //         size_t ownId = neighbours.ownId[i];
    //         size_t neighId = neighbours.neighId[i];

    //         const float val = val_ij[i];
    //         for (int j = 0; j < 3; j++) {
    //             this->f_[ownId][j] += val * kernel.dWdx[i][j];
    //             this->f_[neighId][j] -= val * kernel.dWdx[i][j];
    //         }
    //     }
    // };
};

// NOTE Same as SPHScalarField but different length
template <class T>
class SPHScalarABField: public SPHField<T> {


    public:

        SPHScalarABField(
            std::vector<T>& field,
            const std::string name,
            const std::string type
            ):
            SPHField<T>(name, type, field) {};

    // for(int row = 0; row < 4; row++)
    // {
    //     //calculate_resultant_row(row);
    //     const double* rowA = (const double*)&A[row];
    //     __m256d* pr = (__m256d*)(&C[row]);
    //
    //     *pr = _mm256_mul_pd(_mm256_broadcast_sd(&rowA[0]), matB.row[0]);
    //     for(int i = 1; i < 4; i++)
    //         *pr = _mm256_add_pd(*pr, _mm256_mul_pd(_mm256_broadcast_sd(&rowA[i]),
    //             matB.row[i]));
    // }
};


class SPHIntField: public SPHScalarField<int> {

    public:

        SPHIntField(
            std::vector<int>& field,
            const std::string name=""):
            SPHScalarField<int>(field, name, "int") {};

};

class SPHSizeTField: public SPHScalarField<size_t> {

    public:

        SPHSizeTField(
                std::vector<size_t>& field,
                const std::string name=""):
                SPHScalarField<size_t>(field, name, "long") {};
};

class SPHFloatField: public SPHScalarField<float> {

public:

  SPHFloatField(
                std::vector<float>& field,
                const std::string name = ""
                ) :
    SPHScalarField<float>(field, name, "float") {};

  SPHFloatField(SPHScalarField<float> b)
      : SPHScalarField<float>(b.get_field(), b.get_name(), "float") {};

  // Serial operations

};

template <class T> class SPHComponentField : public SPHField<T> {

  private:
    std::vector<std::string> comp_names_;

  public:
    SPHComponentField(
        std::vector<T> &field,
        const std::string type_name,
        const std::string name,
        std::vector<std::string> comp_names)
      : SPHField<T>(name, type_name, field),
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
        const std::string name = "")
      : SPHComponentField<Vector>(field, "vector", name,  comp_names) {};

  void set(SPHVectorField b) {
    const size_t size = this->f_.size();
    for (size_t ctr = 0; ctr < size; ctr++) {
      this->f_[ctr][0] = b[ctr][0];
      this->f_[ctr][1] = b[ctr][1];
      this->f_[ctr][2] = b[ctr][2];
    }
  }

  SPHVectorField operator*(const float b) const {
    std::vector<Vector> ret(this->f_.size());

    const size_t size = this->f_.size();

    for (size_t i = 0; i < size; i++) {
      ret[i][0] = this->f_[i][0] * b;
      ret[i][1] = this->f_[i][1] * b;
      ret[i][2] = this->f_[i][2] * b;
    }

    return SPHVectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp");
  }

  SPHVectorField operator+(const SPHVectorField b) const {
    std::vector<Vector> ret(this->f_.size());

    const size_t size = this->f_.size();

    for (size_t i = 0; i < size; i++) {
      ret[i][0] = this->f_[i][0] + b.get_field()[i][0];
      ret[i][1] = this->f_[i][1] + b.get_field()[i][1];
      ret[i][2] = this->f_[i][2] + b.get_field()[i][2];
    }

    return SPHVectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp");
  }

};

class SPHPointField : public SPHComponentField<Point> {

  public:
    SPHPointField(std::vector<Point> &field, const std::string name = "")
        : SPHComponentField<Point>(
              field,
              "Point",
              name,
              std::vector<std::string>({"X", "Y", "Z"})) {};

    const Point &operator[](size_t id) const { return this->get_field()[id]; }

  void operator+=(SPHVectorField b) {
    const size_t size = this->f_.size();

    for (size_t i = 0; i < size; i++) {
      this->f_[i] = this->f_[i] + CGALVector(b[i][0], b[i][1], b[i][2]);
    }
  }
};

std::ostream &operator<<(std::ostream &os, SPHFloatField const &m) {

    auto fs = m.get_field();
    os << m.get_name() << std::endl;
    size_t id = 0;
    for (auto f : fs) {
      os << m.get_name() << "Particle " << id << " " << f << "\n" << std::flush;
      id++;
    }
    return os;

}

std::ostream &operator<<(std::ostream &os, SPHVectorField const &m) {

  auto fs = m.get_field();
  os << m.get_name() << std::endl;
  size_t id = 0;
  for (auto f : fs) {os
      << m.get_name() << " Particle " << id << " ("
      << f[0] << " "
      << f[1] << " "
      << f[2] << ")\n" << std::flush;
    id++;
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, SPHPointField const &m) {

  auto fs = m.get_field();
  os << m.get_name() << std::endl;
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

#endif
