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

#endif
