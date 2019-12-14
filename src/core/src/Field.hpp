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
#include "Object.hpp"
#include "cgal/CGALTYPEDEFS.hpp"
#include "SearchCubes.hpp"
#include "Datastructures.hpp"
#include "Vec3.hpp"

#include <iostream>
#include <vector>
#include <math.h>



// A stateful transform that records whether all the std::vector<> terminals
// it has seen are equal to the given size.
struct equal_sizes_impl
{
    template <typename T>
    auto operator() (boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
                     std::vector<T> const & vec)
    {
        auto const expr_size = vec.size();
        if (expr_size != size)
            value = false;
        return 0;
    }

    std::size_t const size;
    bool value;
};

template <typename Expr>
bool equal_sizes (std::size_t size, Expr const & expr)
{
    equal_sizes_impl impl{size, true};
    boost::yap::transform(boost::yap::as_expr(expr), impl);
    return impl.value;
}






// // As assign() above, just using +=.
// template <typename T, typename Expr>
// std::vector<T> & operator+= (std::vector<T> & vec, Expr const & e)
// {
//     decltype(auto) expr = boost::yap::as_expr(e);
//     assert(equal_sizes(vec.size(), expr));
//     for (std::size_t i = 0, size = vec.size(); i < size; ++i) {
//         vec[i] += boost::yap::evaluate(
//             boost::yap::transform(boost::yap::as_expr(expr), take_nth{i}));
//     }
//     return vec;
// }

// As assign() above, just using +=.

template <class T> class Field : public T, public SPHObject {

  private:
    std::vector<std::string> comp_names_;

  protected:
    bool reorder_ = true;

  public:

    Field() {};

    Field(
        T f,
        const std::string name = "",
        std::vector<std::string> comp_names = {})
        : T(f), SPHObject(name, "Field") {
    };

    Field(
        int size,
        typename T::value_type val,
        const std::string name = "",
        std::vector<std::string> comp_names = {})
        : T(size, val), SPHObject(name, "Field") {
    };

    virtual void set_reorder(bool reorder) { reorder_ = reorder; }

    virtual void reorder(const std::vector<size_t> &idx) {
        // // TODO refactor with better error handling system
        // if (*this->size() != idx.size()) {
        //     logger_.warn() << " size mismatch " << name_;
        // }
        // if (reorder_ && *this->size() > 0) {
        //     logger_.info_begin() << __func__ << " reordering " << name_;
        //     reorder_vector(idx, *this);
        //     logger_.info_end();
        // }
    }

    // Setter

    // template <class T, class I>
    // T set(I val);

    // void weighted_sum();

    // void append(std::vector<T> b) {
    //     this->f_.insert(
    //         this->f_.end(),
    //         b.begin(),
    //         b.end()
    //         );
    // };
};

// A derived field to distinguish AB fields from normal fields
//
// A FieldAB is an alias to Field, however, using inheritance
// disallows using Field and FieldAB in binary operators of
// type operator(T a, T b)
template<class T>
class FieldAB: public Field<T> {
    using F = Field<T>;
    using F::F;
};

// A Field wrapper for Field<T>
//
// This field wrapper allows take_nth to dispatch to
// different operator() implementations and select correct
// index
template <class T> class A {
  private:
    Field<T> &f_;

  public:
    A(Field<T> &f) : f_(f) {};

    Field<T> &operator()() { return f_; };

    Field<T> &operator()() const { return f_; };
};

template <class T> class B {

  private:
    Field<T> &f_;

  public:
    B(Field<T> &f) : f_(f) {};

    Field<T> &operator()() const { return f_; };
};

// Free functions

// TODO make it a member function of the fields
template <class T>
void reorder_vector(const std::vector<size_t>& idxs, std::vector<T>& vec);

void write_to_disk(std::string path);

// std::ostream &operator<<(std::ostream &os, Field const &m);

// template<class T>
// T::value_type get_max(const Field<T> &f) {
//     const size_t size = this->f_.size();
//     // TODO implement data type dependent lower bound
//     T ret = std::numeric_limits<T>::min();
//     for (size_t ctr = 0; ctr < size; ctr++) {
//         ret = max(ret, this->f_[ctr]);
//     }
//     return ret;
// }





// Lazy functions

// Define a type trait that identifies std::vectors.
template <typename T>
struct is_field : std::false_type {};

template <typename T, typename A>
struct is_field<Field<std::vector<T, A>>> : std::true_type {};

template <typename T, typename A>
struct is_field<FieldAB<std::vector<T, A>>> : std::true_type {};

template <typename T, typename Alloc>
struct is_field<A<std::vector<T, Alloc>>> : std::true_type {};

// template <>
// struct is_field<Field<std::vector<Vec3>>> : std::true_type {};

template<typename T>
struct zero {};

template<>
struct zero<int> {constexpr static int val = 1;};

template<>
struct zero<float> {constexpr static float val = 0.0;};

template<>
struct zero<Vec3> {constexpr static Vec3 val = {0.0, 0.0, 0.0};};

// Define all the expression-returning numeric operators we need.  Each will
// accept any std::vector<> as any of its arguments, and then any value in the
// remaining argument, if any -- some of the operators below are unary.

// -
BOOST_YAP_USER_UDT_UNARY_OPERATOR(
    negate, boost::yap::expression, is_field);

// *
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    multiplies, boost::yap::expression, is_field);

// /
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    divides, boost::yap::expression, is_field);

// %
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    modulus, boost::yap::expression, is_field);

// +
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    plus, boost::yap::expression, is_field);

// -
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    minus, boost::yap::expression, is_field);

// <
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    less, boost::yap::expression, is_field);

// >
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    greater, boost::yap::expression, is_field);

// <=
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    less_equal, boost::yap::expression, is_field);

// >=
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    greater_equal, boost::yap::expression, is_field);

// ==
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    equal_to, boost::yap::expression, is_field);

// !=
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    not_equal_to, boost::yap::expression, is_field);

// ||
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    logical_or, boost::yap::expression, is_field);


// &&
BOOST_YAP_USER_UDT_ANY_BINARY_OPERATOR(
    logical_and, boost::yap::expression, is_field);


std::ostream &operator<<(std::ostream &os, Vec3 const &f) {
    os << "[" << f[0] << ", " << f[1] << ", " << f[2] << "]";
    return os;
};

template<class T>
std::ostream &operator<<(std::ostream &os, Field<T> const &f) {
    os << f.get_type() << ": " << f.get_name() << "\n[";
    for (auto &v: f) {
        os << "\n " << v;
    }
    os << "\n]\n";
    return os;
};


using FloatField = Field<std::vector<float>>;
using IntField = Field<std::vector<int>>;
using SizeTField = Field<std::vector<size_t>>;

using FloatFieldAB = FieldAB<std::vector<float>>;
using IntFieldAB = FieldAB<std::vector<int>>;
using SizeTFieldAB = FieldAB<std::vector<size_t>>;
using NeighbourFieldAB = FieldAB<std::vector<searchcubes::NeighbourPair>>;
using KernelGradientField = FieldAB<std::vector<VectorPair>>;

using SizeTVectorField = Field<std::vector<std::vector<size_t>>>;

using VectorField = Field<std::vector<Vec3>>;
using PointField = Field<std::vector<Point>>;

// struct to access field elements
//
// This struct implements several operator() versions
// dependent on the passed typed it either uses the
// a-th, b-th, or ab-th index
struct take_nth {
    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        Field<std::vector<T>> const &f) {
        return boost::yap::make_terminal(f[a]);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        FieldAB<std::vector<T>> const &f) {
        return boost::yap::make_terminal(f[ab]);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        A<std::vector<T>> const &f) {
        return boost::yap::make_terminal(f()[a]);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        B<std::vector<T>> const &f) {
        return boost::yap::make_terminal(f()[a]);
    }

    // owner particle index a
    std::size_t a;

    // neighbour particle index b
    std::size_t b;

    // running ab index
    std::size_t ab;
};

// Assigns some expression e to the given vector by evaluating e elementwise,
// to avoid temporaries and allocations.
template <typename T, typename Expr>
std::vector<T> & solve (std::vector<T> & vec, Expr const & e)
{
    decltype(auto) expr = boost::yap::as_expr(e);
    assert(equal_sizes(vec.size(), expr));
    for (std::size_t i = 0, size = vec.size(); i < size; ++i) {
        auto vec_i_expr = boost::yap::transform(boost::yap::as_expr(expr), take_nth{i, 0, 0});
        vec[i] = boost::yap::evaluate(vec_i_expr);
    }
    return vec;
}


// Assigns some expression e to the given vector by evaluating e elementwise,
// to avoid temporaries and allocations.
template <typename T, typename Expr>
std::vector<T> &
sum_AB(std::vector<T> &vec, NeighbourFieldAB &nb, Expr const &e) {
    decltype(auto) expr = boost::yap::as_expr(e);
    // Iterate particle index a
    size_t ab_index = 0;
    for (std::size_t a = 0, size = vec.size(); a < size; ++a) {
        while (a == nb[ab_index].ownId && ab_index < nb.size() ) {
            auto nb_pair = nb[ab_index];
            size_t b = nb_pair.neighId;
            auto vec_ij_expr = boost::yap::transform(
                boost::yap::as_expr(expr), take_nth {a, b, ab_index});
            auto res = boost::yap::evaluate(vec_ij_expr);
            vec[a] += res;
            vec[b] += res;
            ab_index++;
        }
    }
    return vec;
}

// Assigns some expression e to the given vector by evaluating e elementwise,
// to avoid temporaries and allocations.
template <typename T, typename Expr>
std::vector<T> &sum_AB_dW(
    std::vector<T> &vec,
    NeighbourFieldAB &nb,
    KernelGradientField &dW,
    Expr const &e) {
    decltype(auto) expr = boost::yap::as_expr(e);
    // Iterate particle index a
    size_t ab_index = 0;
    for (std::size_t a = 0, size = vec.size(); a < size; ++a) {
        while (a == nb[ab_index].ownId) {
            auto nb_pair = nb[ab_index];
            size_t b = nb_pair.neighId;
            auto vec_ij_expr = boost::yap::transform(
                boost::yap::as_expr(expr), take_nth {a, b, ab_index});
            auto res = boost::yap::evaluate(vec_ij_expr);
            vec[a] += res * dW[ab_index].on;
            vec[b] += res * dW[ab_index].no;
            ab_index++;
        }
    }
    return vec;
}

#endif
