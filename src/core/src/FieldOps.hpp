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

#ifndef PARTIKLER_FIELD_OPS_INCLUDED
#define PARTIKLER_FIELD_OPS_INCLUDED

#include "Field.hpp"
#include "SearchCubes.hpp"
#include <boost/hana/fwd/back.hpp>
#include <boost/yap/print.hpp>
#include <boost/yap/yap.hpp>

#include "Scalar.hpp"
#include <math.h>
#include <queue>
#include <cstdlib>
#include <stdexcept>

template <class ValType> struct Pow_Wrapper {
    Pow_Wrapper() {};
    Pow_Wrapper(ValType i) { inner = i; };
    template <typename T> T operator()(T x) { return std::pow(x, inner); }
    ValType inner;
};

template <class T>
class CyclicVector              {

private:

    std::vector<T> v_;

    // current insert position
    size_t p_;

public:

    CyclicVector(size_t n):
        v_(std::vector<T> (n, zero<T>::val)), p_(0) {
    };

    std::pair<size_t, size_t> div(size_t x, size_t y) {
        size_t rem = x % y;
        size_t quot = x / y;
        return {quot, rem};
    }

    // access operator
    T& operator[](size_t i) {

        auto res = div(p_+i, v_.size());
        // access position should not cross p
        if (res.first > 1 && res.second > p_) {
            throw std::runtime_error("cache error");
        }

        return v_[res.second];
    };

    void add(size_t i, T elem) {

        auto res = div(p_+i, v_.size());
        // access position should not cross p
        if (res.first > 1 && res.second > p_) {
            throw std::runtime_error("cache error");
        }

        v_[res.second] += elem;
    };
    // returns first element and advances it position
    T pop_front() {
        size_t po = p_;
        if (p_ + 1 == v_.size()) {
            p_ = 0;
        } else {
            p_ += 1;
        }
        return v_[po];
    };

};

// template <class OpType, class... Inner> struct Terminal_Generator {
//     Terminal_Generator(Inner... params) {
//         terminal = boost::yap::make_terminal(OpType(params...));
//     }

//     // copy  expressions
//     template <typename T> auto operator()(T x) { return terminal(x); }
//     // but take Fields by reference
//     template <typename U> auto operator()(Field<U> &x) { return terminal(x);
//     }

//     decltype(boost::yap::make_terminal(OpType())) terminal;
// };

// template <class OpType> struct Terminal_Generator_NoArg {
//     Terminal_Generator_NoArg() {
//         terminal = boost::yap::make_terminal(OpType());
//     }
//     // copy  expressions
//     template <typename T> auto operator()(T x) { return terminal(x); }
//     // but take Fields by reference
//     template <typename U> auto operator()(Field<U> &x) { return terminal(x);
//     }

//     decltype(boost::yap::make_terminal(OpType())) terminal;
// };

struct Norm_Wrapper {
    Norm_Wrapper() {};
    template <typename T> Scalar operator()(T x) {
        Scalar ret = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        return ret;
    }
};

struct NormSqr_Wrapper {
    NormSqr_Wrapper() {};
    template <typename T> Scalar operator()(T x) {
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    }
};

template <class ValType> struct Set_Wrapper {
    Set_Wrapper() {};
    Set_Wrapper(ValType i) { inner = i; }
    template <typename T> T operator()(T x) { return inner; }
    ValType inner;
};

struct take_ith {

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>, Scalar const s) {
        return boost::yap::make_terminal(s);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        Field<std::vector<T>> const &f) {
        return boost::yap::make_terminal(f[a]);
    }

    // owner particle index a
    std::size_t a;
};

// struct to access field elements
//
// This struct implements several operator() versions
// dependent on the passed typed it either uses the
// a-th, b-th, or ab-th index or simply makes a terminal
struct take_nth {

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>, Scalar const s) {
        return boost::yap::make_terminal(s);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        Field<std::vector<T>> const &f) {
        return boost::yap::make_terminal(f[a]);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        FieldAB<Field<std::vector<T>>> const &f) {
        return boost::yap::make_terminal(f[ab]);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        A<Field<std::vector<T>>> const &f) {
        return boost::yap::make_terminal(f()[a]);
    }

    template <typename T>
    auto operator()(
        boost::yap::expr_tag<boost::yap::expr_kind::terminal>,
        B<Field<std::vector<T>>> const &f) {
        return boost::yap::make_terminal(f()[b]);
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
std::vector<T> &
solve_impl_res(std::vector<T> &vec, std::vector<bool>& mask, Expr const &e) {
    decltype(auto) expr = boost::yap::as_expr(e);
    std::fill(vec.begin(), vec.end(), zero<T>::val);
    return solve_impl(vec, mask, e);
}

template <typename T, typename Expr>
std::vector<T> &
solve_impl(std::vector<T> &vec, std::vector<bool> &mask, Expr const &e) {
    decltype(auto) expr = boost::yap::as_expr(e);

    size_t thread_begin = 0;
    size_t thread_end = 0;

    for (std::size_t i = 0, size = vec.size(); i < size; ++i) {
        auto vec_i_expr =
            boost::yap::transform(boost::yap::as_expr(expr), take_ith {i});

        // Solve Equation even if this particle is a boundary particle
        // since evaluate sets the neighbour particle value as side effect

        auto res = boost::yap::evaluate(vec_i_expr);

        // check if i is a boundary
        // NOTE
        // if fixed value do nothing since value should be assigned on field
        // construction

        if ( mask[i] ) {
            vec[i] = res;
        }

    }
    return vec;
}

// appends the values of b to a inplace
template <class T> void field_append(T &a, T &b) {
    a.insert(a.end(), b.begin(), b.end());
}

template <class Inner> struct Sum_AB_sym {

    Sum_AB_sym(Inner &state, NeighbourFieldAB &nb)
        : vec_(state),
          cache_(state.size()),
          nb_(nb)
    {};

    template <class Expr> typename Inner::value_type operator()(Expr expr) {

        // load orig a value
        // TODO pass a flag if needed use zero otherwise, this would make
        // resetting unneeded

        ret_ = vec_[a];

        // if a cached operation is present
        // apply and pop it
        ret_ += cache_.pop_front();

        while (ab < nb_.size() && a == nb_[ab].ownId) {
            auto nb_pair = nb_[ab];
            size_t b = nb_pair.neighId;
            size_t dist = b-a-1;

            auto vec_ij_expr = boost::yap::transform(
                boost::yap::as_expr(expr), take_nth {a, b, ab});
            auto res = boost::yap::evaluate(vec_ij_expr);
            ret_ += res;

            cache_.add(dist, res);
            ab++;
        }
        a++;
        return ret_;
    }

    size_t a = 0;
    size_t ab = 0;
    typename Inner::value_type ret_= zero<typename Inner::value_type>::val;

    Inner &vec_;
    CyclicVector <typename Inner::value_type> cache_;

    NeighbourFieldAB &nb_;
};

template <class Inner> struct Sum_AB_asym {
    Sum_AB_asym(Inner &state, NeighbourFieldAB &nb)
        : vec_(state),
          cache_(state.size()),
          nb_(nb)
    {};

    template <class Expr> typename Inner::value_type operator()(Expr expr) {

        // load orig a value
        // TODO pass a flag if needed use zero otherwise, this would make
        // resetting unneeded

        ret_ = vec_[a];

        // if a cached operation is present
        // apply and pop it
        ret_ -= cache_.pop_front();

        while (ab < nb_.size() && a == nb_[ab].ownId) {
            auto nb_pair = nb_[ab];
            size_t b = nb_pair.neighId;
            size_t dist = b-a-1;

            auto vec_ij_expr = boost::yap::transform(
                boost::yap::as_expr(expr), take_nth {a, b, ab});
            auto res = boost::yap::evaluate(vec_ij_expr);
            ret_ += res;

            cache_.add(dist, res);
            ab++;
        }
        a++;
        return ret_;
    }

    size_t a = 0;
    size_t ab = 0;
    typename Inner::value_type ret_= zero<typename Inner::value_type>::val;

    Inner &vec_;
    CyclicVector <typename Inner::value_type> cache_;

    NeighbourFieldAB &nb_;
};

template <class Inner> struct Sum_AB_dW_sym {

    Sum_AB_dW_sym(
        Inner &state,
        NeighbourFieldAB &nb,
        KernelGradientField &dW,
        KernelGradientField &dWn)
        : vec_(state), cache_(state.size()), nb_(nb), dW_(dW), dWn_(dWn) {};

    template <class Expr> typename Inner::value_type operator()(Expr expr) {

        // load orig a value
        // TODO pass a flag if needed use zero otherwise, this would make
        // resetting unneeded

        ret_ = vec_[a];

        // if a cached operation is present
        // apply and pop it
        ret_ -= cache_.pop_front();

        while (ab < nb_.size() && a == nb_[ab].ownId) {
            auto nb_pair = nb_[ab];
            size_t b = nb_pair.neighId;
            size_t dist = b-a-1;
            auto vec_ij_expr = boost::yap::transform(
                boost::yap::as_expr(expr), take_nth {a, b, ab});
            auto res = boost::yap::evaluate(vec_ij_expr);
            ret_ += res * dW_[ab];

            // rather use a vector since to keep size constant
            cache_.add(dist, res*dWn_[ab]);
            ab++;
        }
        a++;
        return ret_;
    }

    size_t a = 0;
    size_t ab = 0;
    typename Inner::value_type ret_= zero<typename Inner::value_type>::val;

    Inner &vec_;
    CyclicVector <typename Inner::value_type> cache_;
    NeighbourFieldAB &nb_;
    KernelGradientField &dW_;
    KernelGradientField &dWn_;
};

template <class Inner> struct Sum_AB_dW_asym {

    Sum_AB_dW_asym(
        Inner &state,
        NeighbourFieldAB &nb,
        KernelGradientField &dW,
        KernelGradientField &dWn)
        : vec_(state), cache_(state.size()), nb_(nb), dW_(dW), dWn_(dWn) {};

    template <class Expr> typename Inner::value_type operator()(Expr expr) {
        // load orig a value
        // TODO pass a flag if needed use zero otherwise, this would make
        // resetting unneeded

        ret_ = vec_[a];

        // if a cached operation is present
        // apply and pop it
        ret_ += cache_.pop_front();

        while (ab < nb_.size() && a == nb_[ab].ownId) {
            auto nb_pair = nb_[ab];
            size_t b = nb_pair.neighId;
            size_t dist = b-a-1;
            auto vec_ij_expr = boost::yap::transform(
                boost::yap::as_expr(expr), take_nth {a, b, ab});
            auto res = boost::yap::evaluate(vec_ij_expr);
            ret_ += res * dW_[ab];

            // rather use a vector since to keep size constant
            cache_.add(dist, res*dWn_[ab]);
            ab++;
        }
        a++;
        return ret_;
    }

    size_t a = 0;
    size_t ab = 0;
    typename Inner::value_type ret_= zero<typename Inner::value_type>::val;


    Inner &vec_;
    CyclicVector <typename Inner::value_type> cache_;
    NeighbourFieldAB &nb_;
    KernelGradientField &dW_;
    KernelGradientField &dWn_;
};

template <typename T, typename Expr>
std::vector<T> &
solve_inner_impl(NeighbourFieldAB &nb, std::vector<T> &vec, Expr const &e) {
    decltype(auto) expr = boost::yap::as_expr(e);
    for (std::size_t i = 0, size = vec.size(); i < size; ++i) {
        auto nb_pair = nb[i];
        size_t b = nb_pair.neighId;
        size_t a = nb_pair.ownId;
        auto vec_ij_expr = boost::yap::transform(
            boost::yap::as_expr(expr), take_nth {a, b, 0});
        auto res = boost::yap::evaluate(vec_ij_expr);

        vec[i] = res;
    }
    return vec;
};

template <class Inner> struct Ddt {

    // // Default constructor to allow decltype(make_terminal(Ddt()))
    // construction Ddt() {};

    Ddt(Scalar dt, Inner &vec)
        : dt_(dt), vec_({}) {
        // store old state
        vec_.reserve(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            vec_.push_back(vec[i]);
        }
    };

    template <class T> T operator()(T v) {

        T ret = vec_[a_] + v * dt_;
        a_++;
        return ret;
    }

    size_t a_ = 0;

    Scalar dt_;

    std::vector<typename Inner::value_type> vec_;
};

template <class Field, class El>
void clamp_field_in_range(Field &f, El lo, El hi);

#endif
