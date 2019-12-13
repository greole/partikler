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
    {return std::pow(x, inner);};
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
    {return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);};
};

template<class... ValType>
using Norm = Terminal_Generator<Norm_Wrapper, ValType...>;

template<class ValType>
struct Set_Wrapper {
    Set_Wrapper() {};
    Set_Wrapper(ValType i) { inner = i; }
    template <typename T> T operator()(T x)
    {
        return inner;};
    ValType inner;
};

template<class... ValType>
using Set = Terminal_Generator<Set_Wrapper<ValType...>, ValType...>;

// struct Sum_Wrapper {
//     Sum_Wrapper(SizeTVectorField & nb, FloatField& f):nb_(nb),f_(f)  {
//     };
//     template <typename T> float operator()(T x) {
//         float ret = 0;
//         for (auto i : nb_[x]) {
//             ret += f_[i];
//         }

//         return ret;
//     };

//     SizeTVectorField& nb_;
//     FloatField& f_;
// };

// struct Weighted_Sum_Wrapper {

//     Sum_Wrapper(SizeTVectorField& nb, FloatField& f):nb_(nb), f_(f)  {};

//     template <typename T> float operator()(T x) {
//         float ret = 0;
//         for (auto i : nb_[x]) {
//             ret += f_[i];
//         }

//         return ret;
//     };

//     SizeTVectorField& nb_;
//     FloatVectorField& W_;
//     FloatField& f_;
// };



// template <class T>
// Field<std::vector<T>>
// sum_AB(Field<std::vector<T>> f, NeighbourFieldAB nbs, FloatFieldAB W) {

//     Field<std::vector<T>> ret(f.size(), zero<T>::val);

//     for (size_t i = 0; i < W.size(); i++) {
//         ret[nbs[i].ownId] += f*W[i];
//         ret[nbs[i].neighId] += f*W[i];
//     }
//     return ret;
// };

// template <class T>
// struct sum_AB {

//     sum_AB(size_t n, NeighbourFieldAB &nbs, FloatFieldAB &W)
//         : n_(n),nbs_(nbs), W_(W) {};

//     auto operator()() {
//         Field<std::vector<T>> ret(n_, zero<T>::val);

//         for (size_t i = 0; i < W_.size(); i++) {
//             ret[nbs_[i].ownId] += W_[i];
//             ret[nbs_[i].neighId] += W_[i];
//         }

//         return ret;
//     }

//     size_t n_;
//     NeighbourFieldAB &nbs_;
//     FloatField &W_;
// };

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
};


    // template <class T>
    // Field<T> sum_AB(
    //     Field<T> f,
    //     FieldAB<searchcubes::NeighbourPair> nbs,
    //     FieldAB<std::vector<float>> W) {

    //     Field<T> ret ();

    //     for (auto nb:nbs) {
    //         ret[nb.first] += ;
    //         ret[nb.second] +=
    //     }

    // };

    // template <class T>
    // Field<T> sum_AB(
    //     Field<T> f,
    //     FieldAB<searchcubes::NeighbourPair> nb,
    //     FieldAB<VectorPair> dWdx) {

    // };

#endif
