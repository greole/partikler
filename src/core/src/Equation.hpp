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
#include "SearchCubes.hpp"    // for NeighbourFieldAB
#include "Time.hpp"
#include "yaml-cpp/node/impl.h" // for Node::~Node, Node::Node, Node::as
#include "yaml-cpp/node/node.h" // for Node
#include "yaml-cpp/yaml.h"

// TODO TransportedQuantity
// a class providing field to Equation depedency
// has an update method that calls the dependend transport eqn

// TODO TransportEqn Base class
// - Access previous timestep result if activated
// - has an update method
// - allows other TransportEqn as dependencies
// - Computes a d(...)

// class EmptyFieldEquation {

// };

template <class FieldGradientType> class FieldGradientEquation {

  protected:
    FieldGradientType &df_;

  public:
    typedef FieldGradientType value_type;

    template <class FieldType>
    FieldGradientEquation(ObjectRegistry &objReg, FieldType &f)
        : df_(objReg.create_field<FieldGradientType>(
              "d" + f.get_name(),
              zero<typename FieldGradientType::value_type>::val,
              {"dx" + f.get_name(),
               "dy" + f.get_name(),
               "dz" + f.get_name()})) {}
};

template <class FieldType> class FieldValueEquation {

  protected:
    FieldType &f_;

  public:
    typedef FieldType value_type;

    FieldValueEquation(ObjectRegistry &objReg, FieldType &f) : f_(f) {};
};

template <class FieldValueType, class FieldGradientType>
class FieldEquationBase : public Model,
                          public FieldValueType,
                          public FieldGradientType {

  protected:
    TimeGraph &time_;

    NeighbourFieldAB &np_;

    ScalarField &W_;

    // store previous results
    // std::map<int, T> prev_;

    // the current iteraton
    int iteration_;

  public:
    KernelGradientType kernelGradientType_;

    template <class FieldType>
    FieldEquationBase(
        const std::string &field_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        FieldType &f)
        : Model(field_name, parameter, objReg), FieldValueType(objReg, f),
          FieldGradientType(objReg, f),
          time_(objReg.get_object<TimeGraph>("TimeGraph")),
          np_(objReg.get_object<NeighbourFieldAB>("neighbour_pairs")),
          W_(objReg.get_object<ScalarField>("KernelW")),
          // TODO refactor kernelGradientType to FieldGradientType
          kernelGradientType_(
              (read_or_default_coeff_impl<std::string>(
                   parameter, "KernelType", "Symmetric") == "Symmetric")
                  ? KernelGradientType::Symmetric
                  : KernelGradientType::NonSymmetric) {}

    // get result for iteration i
    // if result is not cached solve gets executed

    ScalarField &W() { return W_; };

    NeighbourFieldAB &N() { return np_; };

    template <class Expr>
    typename FieldGradientType::value_type &
    sum_AB_dW(ScalarField &rho, Expr const &e) {
        // dispatch on KernelGradientType
        auto &dW =
            get_objReg().get_object<KernelGradientField>("KerneldWdx");
        decltype(auto) expr = boost::yap::as_expr(e);
        sum_AB_dW_res_impl_rho(rho, this->df_, np_, dW, expr);

        return this->df_;
    }

    template <class Expr>
    typename FieldGradientType::value_type &sum_AB_dW(Expr const &e) {
        // dispatch on KernelGradientType
        if (kernelGradientType_ == KernelGradientType::Symmetric) {
            auto &dW =
                get_objReg().get_object<KernelGradientField>("KerneldWdx");
            decltype(auto) expr = boost::yap::as_expr(e);
            sum_AB_dW_res_impl(this->df_, np_, dW, expr);
        } else {
            auto &dW = get_objReg().get_object<DoubleKernelGradientField>(
                "KerneldWdx");
            decltype(auto) expr = boost::yap::as_expr(e);
            sum_AB_dW_res_impl(this->df_, np_, dW, expr);
        }

        return this->df_;
    }

    template <class RHS> typename FieldValueType::value_type &sum_AB(RHS rhs) {
        sum_AB_res_impl(this->f_, np_, rhs * W_);
        return this->f_;
    }

    typename FieldValueType::value_type &m_sum_AB(float particle_mass) {
        sum_AB_res_impl(particle_mass, this->f_, np_, W_);
        return this->f_;
    };

    // VectorField& ddx() {
    // };

    template <class RHS> void ddt(RHS rhs) {

        // integrate RK4 -> k1 = ddt(t, u^o).get(), k2=ddt(t+h/2, u^o+h/2*k1)
        // ...
    }

    // template<class T>
    // virtual T expression();

    // solve();

    // get the result for the given iteration
    //
    // get the result for the given iteration
    // if iteration is larger then cached iteration
    // equation is solved
    typename FieldValueType::value_type &get(int iteration) {
        if (iteration == iteration_) {
            std::cout << get_name()
                      << "DEBUG using cached solution for iteration "
                      << iteration << std::endl;
            return this->f_;
        } else {
            this->execute();
            return this->f_;
        }
    }

    typename FieldGradientType::value_type &get_dx(int iteration) {
        if (iteration == iteration_) {
            std::cout << get_name()
                      << "DEBUG using cached solution for iteration "
                      << iteration << std::endl;
            return this->df_;
        } else {
            this->execute();
            return this->df_;
        }
    }
};

using ScalarFieldEquation = FieldEquationBase<
    FieldValueEquation<ScalarField>,
    FieldGradientEquation<VectorField>>;

using VectorFieldEquation = FieldEquationBase<
    FieldValueEquation<VectorField>,
    FieldGradientEquation<VectorField>>;

// TODO implement this
// using ABVectorFieldEquation =
//         FieldEquation<VectorFieldAB, VectorFieldAB>;

#endif
