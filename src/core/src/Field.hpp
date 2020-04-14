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
#ifndef PARTIKLER_FIELD_INCLUDED
#define PARTIKLER_FIELD_INCLUDED

#include <boost/yap/yap.hpp>
#include <stddef.h>    // for size_t
#include <string>      // for string, operator<<
#include <type_traits> // for true_type, false_type
#include <vector>      // for vector, allocator

#include "Object.hpp"            // for SPHObject, GetFieldType, FloatF...
#include "Scalar.hpp"            // for Vec3, VectorPair (ptr only)
#include "Vec3.hpp"              // for Vec3, VectorPair (ptr only)

#define ab(field) (field.a() - field.b())

// Dynamically dispatches func based on its kind
#define DISPATCH(obj, func, type_enum, ...)                                    \
    switch (type_enum) {                                                       \
    case IntFieldType:                                                         \
        func<IntField>(dynamic_cast<IntField &>(*obj->get()), __VA_ARGS__);    \
        break;                                                                 \
    case SizeTFieldType:                                                       \
        func<SizeTField>(                                                      \
            dynamic_cast<SizeTField &>(*obj->get()), __VA_ARGS__);             \
        break;                                                                 \
    case ScalarFieldType:                                                      \
        func<ScalarField>(                                                     \
            dynamic_cast<ScalarField &>(*obj->get()), __VA_ARGS__);            \
        break;                                                                 \
    case VectorFieldType:                                                      \
        func<VectorField>(                                                     \
            dynamic_cast<VectorField &>(*obj->get()), __VA_ARGS__);            \
        break;                                                                 \
    default:                                                                   \
        continue;                                                              \
    }

enum ReadOptions { MUST_READ, NO_READ };

enum WriteOptions { NO_WRITE, WRITE };

struct IOOptions {
    ReadOptions r;
    WriteOptions w;
};

// A derived field to distinguish AB fields from normal fields
//
// A FieldAB is an alias to Field, however, using inheritance
// disallows using Field and FieldAB in binary operators of
// type operator(T a, T b)
template <class T> class FieldAB : public T {
    using F = T;
    using F::F;
};

// A Field wrapper for Field<T>
//
// This field wrapper allows take_nth to dispatch to
// different operator() implementations and select correct
// index
template <class T> class A {
  private:
    T &f_;

  public:
    A(T &f) : f_(f) {};

    T &operator()() { return f_; };

    T &operator()() const { return f_; };
};

template <class T> class B {

  private:
    T &f_;

  public:
    B(T &f) : f_(f) {};

    T &operator()() const { return f_; };
};

template <class T> class Field : public T, public SPHObject {

  private:
    std::vector<std::string> comp_names_;

    IOOptions iooptions_;

  protected:
    bool reorder_ = true;

  public:
    Field() {};

    Field(
        T f,
        const std::string name = "",
        std::vector<std::string> comp_names = {})
        : T(f), SPHObject(name, GetFieldType<T>::value),
          comp_names_(comp_names) {};

    Field(
        size_t size,
        typename T::value_type val,
        const std::string name = "",
        const std::vector<std::string> comp_names = {})
        : T(size, val), SPHObject(name, GetFieldType<T>::value),
          comp_names_(comp_names) {};

    void set_io_options(IOOptions ioo) { iooptions_ = ioo; };

    IOOptions get_io_options() { return iooptions_; };

    void operator=(Field<T> &b) { T::operator=(b); }

    void operator=(T &b) { T::operator=(b); }

    void operator=(Field<T> &&b) { T::operator=(std::move(b)); }

    void set_reorder(bool reorder) { reorder_ = reorder; }

    A<Field<T>> a() { return A(*this); }

    B<Field<T>> b() { return B(*this); }

    // reoder the vector by the idx vector
    void reorder(const std::vector<size_t> &idx) {
        // // TODO refactor with better error handling system
        if (T::size() != idx.size()) {
            logger_.warn() << " size mismatch " << name_;
        }
        if (reorder_ && T::size() > 0) {
            logger_.info_begin()
                // << __PRETTY_FUNCTION__
                << " reordering " << name_;
            reorder_vector(*this, idx);
            logger_.info_end();
        }
    }

    std::vector<std::string> get_comp_names() const { return comp_names_; };
};

// free function re-oders the vector by the idx vector
// if it is anything ie generic or models do nothing
// partial template specialisation takes care of the fields
template <class T>
void reorder_vector(Field<T> &vec, const std::vector<size_t> &idxs) {
    std::cout << " reorder vector begin" << std::endl;
    T tmp(vec.size());
    for (size_t i = 0; i < idxs.size(); i++) {
        tmp[idxs[i]] = vec[i];
    }

    vec = tmp;
    std::cout << " done" << std::endl;
}

template <class T>
void reorder_vector(
    std::vector<typename T::value_type> &vec, const std::vector<size_t> &idxs) {
    std::vector<typename T::value_type> tmp(vec.size());
    std::cout << " reorder vector II begin" << std::endl;

    for (size_t i = 0; i < idxs.size(); i++) {
        tmp[idxs[i]] = vec[i];
    }

    vec = tmp;
    std::cout << " done" << std::endl;
}

// Lazy functions

// Define a type trait that identifies std::vectors.
template <typename T> struct is_field : std::false_type {};

template <typename T, typename Alloc>
struct is_field<Field<std::vector<T, Alloc>>> : std::true_type {};

template <typename T, typename Alloc>
struct is_field<FieldAB<Field<std::vector<T, Alloc>>>> : std::true_type {};

template <typename T, typename Alloc>
struct is_field<A<Field<std::vector<T, Alloc>>>> : std::true_type {};

template <typename T, typename Alloc>
struct is_field<B<Field<std::vector<T, Alloc>>>> : std::true_type {};

// template <>
// struct is_field<Field<std::vector<Vec3>>> : std::true_type {};

template <typename T> struct zero {};

template <> struct zero<int> { constexpr static int val = 1; };

template <> struct zero<Scalar> { constexpr static Scalar val = 0.0; };

template <> struct zero<Vec3> {
    constexpr static Vec3 val = {{0.0, 0.0, 0.0}};
};

// Define all the expression-returning numeric operators we need.  Each will
// accept any std::vector<> as any of its arguments, and then any value in the
// remaining argument, if any -- some of the operators below are unary.

// -
BOOST_YAP_USER_UDT_UNARY_OPERATOR(negate, boost::yap::expression, is_field)

// *
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    multiplies, boost::yap::expression, is_field)

// /
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    divides, boost::yap::expression, is_field)

// %
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    modulus, boost::yap::expression, is_field)

// +
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(plus, boost::yap::expression, is_field)

// -
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(minus, boost::yap::expression, is_field)

// <
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(less, boost::yap::expression, is_field)

// >
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    greater, boost::yap::expression, is_field)

// <=
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    less_equal, boost::yap::expression, is_field)

// >=
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    greater_equal, boost::yap::expression, is_field)

// ==
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    equal_to, boost::yap::expression, is_field)

// !=
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    not_equal_to, boost::yap::expression, is_field)

// ||
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    logical_or, boost::yap::expression, is_field)

// &&
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    logical_and, boost::yap::expression, is_field)

#include "SearchCubes.hpp"

template <class T>
std::ostream &operator<<(std::ostream &os, Field<T> const &f) {
    os << f.get_type() << ": " << f.get_name() << "\n[";
    size_t it = 0;
    for (auto &v : f) {
        os << "\n " << it << " " << v;
        it++;
    }
    os << "\n]\n";
    return os;
}

using FloatField = Field<std::vector<float>>;
using DoubleField = Field<std::vector<double>>;
using ScalarField = Field<std::vector<Scalar>>;
using IntField = Field<std::vector<int>>;
using SizeTField = Field<std::vector<size_t>>;
using VectorField = Field<std::vector<Vec3>>;
using DoubleVectorField = Field<std::vector<VectorPair>>;

using ScalarFieldAB = FieldAB<ScalarField>;
using VectorFieldAB = FieldAB<Field<std::vector<Vec3>>>;
using IntFieldAB = FieldAB<IntField>;
using SizeTFieldAB = FieldAB<SizeTField>;
using KernelGradientField = FieldAB<VectorField>;
using DoubleKernelGradientField = FieldAB<DoubleVectorField>;

using SizeTVectorField = Field<std::vector<std::vector<size_t>>>;

VectorField &operator+=(VectorField &a, VectorField &b);
VectorField &operator-=(VectorField &a, VectorField &b);

// Template meta function to get SPHObjectType from std::vector<T>
template <enum SPHObjectType T> struct GetField {};

template <> struct GetField<ScalarFieldType> { using type = ScalarField; };

#endif
