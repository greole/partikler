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

#ifndef PARTIKLER_EQUATION_INCLUDED
#define PARTIKLER_EQUATION_INCLUDED

#include <algorithm>                // for max
#include <boost/hana/at.hpp>        // for at_t::operator()
#include <boost/yap/expression.hpp> // for as_expr, make_terminal
#include <iostream>                 // for operator<<, basic_ostream, endl
#include <map>                      // for map<>::iterator, operator!=, ope...
#include <memory>                   // for shared_ptr, allocator, __shared_...
#include <stdio.h>
#include <string>  // for string, operator<<, operator+
#include <utility> // for move, pair, make_pair
#include <vector>  // for vector

#include "Field.hpp"    // for ScalarField, KernelGradientField
#include "FieldOps.hpp" // for sum_AB_impl
#include "Logger.hpp"   // for MSG, Logger
#include "Models.hpp"
#include "Object.hpp"         // for SPHObject, ModelType
#include "ObjectRegistry.hpp" // for ObjectRegistry
#include "Scalar.hpp"
#include "SearchCubes.hpp" // for NeighbourFieldAB
#include "Time.hpp"
#include "yaml-cpp/yaml.h"

// TODO TransportedQuantity
// a class providing field to Equation depedency
// has an update method that calls the dependend transport eqn

// TODO TransportEqn Base class
// - Access previous timestep result if activated
// - has an update method
// - allows other TransportEqn as dependencies
// - Computes a d(...)

template <class FieldType>
class FieldEquationBase : public Model {

  protected:

    FieldType &f_;

    TimeGraph &time_;

    NeighbourFieldAB &np_;

    ObjectRegistry &objReg_;

    // store previous results
    // std::map<int, T> prev_;

    // the current iteraton
    int iteration_;

    Scalar h_;

    Scalar maxDt_  = std::numeric_limits<Scalar>::max();

    Ddt<FieldType> ddt_;

    decltype(boost::yap::make_terminal(Ddt<FieldType>())) ddt_terminal_;

  public:

    FieldEquationBase(
        const std::string &field_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        FieldType& f)
        : Model(field_name, parameter, objReg),
          f_(f),
          time_(objReg.get_object<TimeGraph>("TimeGraph")),
          np_(objReg.get_object<NeighbourFieldAB>("neighbour_pairs")),
          h_(objReg.get_object<Generic<Scalar>>("h")()),
          objReg_(objReg),
          ddt_(Ddt(time_.get_deltaT(),  f_))
           {

          ddt_terminal_ = boost::yap::make_terminal(ddt_);
}

    std::pair<bool, Scalar> timestep_limit() { return {false, 0}; }

    template <class RHS> auto ddt(RHS rhs) {
        return ddt_terminal_(rhs);
    }

    template<class Expr>
    void solve(Expr e) {

        this->log().info_begin() << "Solving: " << f_.name();
        decltype(auto) expr = boost::yap::as_expr(e);
        solve_impl(f_, expr);
        this->log().info_end();
    }

    FieldType &get(int iteration) {
        if (iteration == iteration_) {
            std::cout << get_name()
                      << "DEBUG using cached solution for iteration "
                      << iteration << std::endl;
            return f_;
        } else {
            this->execute();
            return f_;
        }
    }

    // template<class T>
    // virtual T expression();

    // solve();

    // get the result for the given iteration
    //
    // get the result for the given iteration
    // if iteration is larger then cached iteration
    // equation is solved
};

// template <class FieldType>
// class TimeDerivativeEquation : public FieldEquationBase<FieldType> {

// protected:

//     FieldType& Intf_;

// public:

//     TimeDerivativeEquation(
//         const std::string &eqn_base_name,
//         YAML::Node parameter,
//         ObjectRegistry &objReg,
//         FieldType &f):
//          FieldEquationBase<FieldType>(
//             eqn_base_name,
//             parameter,
//             objReg,
//             objReg.create_field<FieldType>(
//                 "d" + f.get_name(),
//                 zero<typename FieldType::value_type>::val,
//                 {"d_" + f.get_name() + "dx",
//                  "d_" + f.get_name() + "dy",
//                  "d_" + f.get_name() + "dz"})),
//          Intf_(f) {}


//     void IntDt() {
//         Scalar dt =  this->time_.get_deltaT();
//         // clang-format off
//         this->log().info_begin()
//             << "Time integration " << Intf_.get_name()
//             << " deltaT " << dt;
//         // clang-format on
//         for (size_t i = 0; i < this->f_.size(); i++) {
//             Intf_[i] += this->f_[i] * dt;
//         }
//         this->log().info_end();
//     }

//     // template <class Expr>
//     // FieldGradientType &sum_AB_dW(Expr const &e) {
//     //     // dispatch on KernelGradientType
//     //     auto &dW = this->get_objReg().template get_object<KernelGradientField>(
//     //         "KerneldWdx");
//     //     decltype(auto) expr = boost::yap::as_expr(e);
//     //     sum_AB_dW_res_impl(rho, this->f_, this->np_, dW, expr);

//     //     return this->f_;
//     // }

// };


template <class FieldGradientType>
class FieldGradientEquation : public FieldEquationBase<FieldGradientType> {

    // TODO use two separate fields for gradients of the Kernel
    // in the symmetric case the neighbour owner gradient field is just a reference
    // other wise the two fields are different fields

  public:
    KernelGradientType kernelGradientType_;

    typedef FieldGradientType value_type;

    template <class FieldType>
    FieldGradientEquation(
        const std::string &eqn_base_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        FieldType &f)
        : FieldEquationBase<FieldGradientType>(
              eqn_base_name,
              parameter,
              objReg,
              objReg.create_field<FieldGradientType>(
                  "d" + f.get_name(),
                  zero<typename FieldGradientType::value_type>::val,
                  {"d_" + f.get_name() + "_dx",
                   "d_" + f.get_name() + "_dy",
                   "d_" + f.get_name() + "_dz"})),
          kernelGradientType_(
              (read_or_default_coeff_impl<std::string>(
                   parameter, "KernelType", "Symmetric") == "Symmetric")
                  ? KernelGradientType::Symmetric
                  : KernelGradientType::NonSymmetric) {}

    template <class Expr>
    FieldGradientType &sum_AB_dW(ScalarField &rho, Expr const &e) {
        // dispatch on KernelGradientType
        auto &dW = this->get_objReg().template get_object<KernelGradientField>(
            "KerneldWdx");
        decltype(auto) expr = boost::yap::as_expr(e);
        sum_AB_dW_res_impl_rho(rho, this->f_, this->np_, dW, expr);

        return this->f_;
    }

    template <class Expr> FieldGradientType &sum_AB_dW(Expr const &e) {
        // dispatch on KernelGradientType
        if (kernelGradientType_ == KernelGradientType::Symmetric) {
            auto &dW = this->get_objReg().template get_object<KernelGradientField>(
                "KerneldWdx");
            decltype(auto) expr = boost::yap::as_expr(e);
            sum_AB_dW_res_impl(this->f_, this->np_, dW, expr);
        } else {
            auto &dW =
                this->get_objReg().template get_object<DoubleKernelGradientField>(
                    "KerneldWdx");
            decltype(auto) expr = boost::yap::as_expr(e);
            sum_AB_dW_res_impl(this->f_, this->np_, dW, expr);
        }

        return this->f_;
    }

    template <class RHS> void solve(RHS rhs, bool reset=false) {
        if (reset) solve_impl_res(this->f_, rhs);
        else solve_impl(this->f_, rhs);
    }

};

template<class FieldType, template<typename> class SumType>
class FieldValueEquation : public FieldEquationBase<FieldType> {

  protected:

    ScalarField &W_;

    SumType<FieldType> sum_AB_s;

  public:
    typedef FieldType value_type;

    FieldValueEquation(
        const std::string &eqn_base_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        FieldType &f)
        : FieldEquationBase<FieldType>(eqn_base_name, parameter, objReg, f),
          W_(objReg.get_object<ScalarField>("KernelW")),
          sum_AB_s(SumType<FieldType>(this->f_, this->np_)) {};

    template <class RHS> void solve(RHS rhs, bool reset=false) {
        if (reset) solve_impl_res(this->f_, rhs);
        else solve_impl(this->f_, rhs);
        // TODO find a better solution
        // Reset sum_AB indizes
        sum_AB_s.a = 0;
        sum_AB_s.ab = 0;
        this->ddt_.a_ = 0;
    }


};


using ScalarFieldEquation = FieldValueEquation<ScalarField, Sum_AB_sym>;
using VectorFieldEquation = FieldValueEquation<VectorField, Sum_AB_sym>;

using ScalarGradientEquation = FieldGradientEquation<VectorField>;
// TODO clarify why this is not a tensor field
using VectorGradientEquation = FieldGradientEquation<VectorField>;


#endif
