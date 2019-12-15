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

#include <boost/yap/yap.hpp>
#include "Field.hpp"
#include "SearchCubes.hpp"


template<class ValType>
struct Pow_Wrapper {
    Pow_Wrapper() {};
    Pow_Wrapper(ValType i) { inner = i; };
    template <typename T> T operator()(T x)
    {return std::pow(x, inner);}
    ValType inner;
};

template <class OpType, class... Inner>
struct Terminal_Generator {
    Terminal_Generator(Inner... params) {
        terminal = boost::yap::make_terminal(OpType (params...));
    }
    template <typename T> auto operator()(T x) { return terminal(x); }

    decltype(boost::yap::make_terminal(OpType())) terminal;
};

template<class... ValType>
using Pow = Terminal_Generator<Pow_Wrapper<ValType...>, ValType...>;


struct Norm_Wrapper {
    Norm_Wrapper() {};
    template <typename T> float operator()(T x)
    {return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);}
};

template<class... ValType>
using Norm = Terminal_Generator<Norm_Wrapper, ValType...>;

template<class ValType>
struct Set_Wrapper {
    Set_Wrapper() {};
    Set_Wrapper(ValType i) { inner = i; }
    template <typename T> T operator()(T x)
    {
        return inner;}
    ValType inner;
};

template<class... ValType>
using Set = Terminal_Generator<Set_Wrapper<ValType...>, ValType...>;


template <boost::yap::expr_kind Kind, typename eqn>
struct equation
{
    static const boost::yap::expr_kind kind = Kind;

    eqn elements;
};

template<class T, class I>
T solve(I terminal) {
    auto expr = boost::yap::make_expression<
        equation,
        boost::yap::expr_kind::call
        >(terminal);

    std::cout << "solving " << std::endl;
    auto res =  boost::yap::evaluate(expr);
    std::cout << "end " << std::endl;

    return res;
}

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
