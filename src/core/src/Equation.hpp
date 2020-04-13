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
#include <boost/yap/print.hpp>      // for as_expr, make_terminal
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

template <class FieldType> class FieldEquationBase : public Model {

  protected:
    FieldType &f_;

    // old time step value
    FieldType fo_;

    IntField &id_;

    TimeGraph &time_;

    NeighbourFieldAB &np_;

    Scalar h_;

    ObjectRegistry &objReg_;

    bool estimate_;

    // store previous results
    // std::map<int, T> prev_;

    Scalar maxDt_ = std::numeric_limits<Scalar>::max();

    int stored_iteration_;

    int iteration_;

    // indicate if strict bounds are enforced
    bool bounded_ {false};
    std::pair<Scalar, Scalar> bounds_;

  public:
    FieldEquationBase(
        const std::string &field_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        FieldType &f)
        : Model(field_name, parameter, objReg), f_(f),
          fo_(f.size(), zero<typename FieldType::value_type>::val),
          id_(objReg.get_object<IntField>("id")),
          time_(objReg.get_object<TimeGraph>("TimeGraph")),
          np_(objReg.get_object<NeighbourFieldAB>("neighbour_pairs")),
          h_(objReg.get_object<Generic<Scalar>>("h")()), objReg_(objReg),
          estimate_(false) {}

    void store_old_value() {
        for (size_t i = 0; i < fo_.size(); i++) {
            fo_[i] = f_[i];
        }
    }

    std::pair<bool, Scalar> timestep_limit() { return {false, 0}; }

    Ddt<FieldType> ddt() { return Ddt(time_.get_deltaT(), fo_, this->id_); }

    template <class Expr> void solve(Expr e) {
        this->log().info_begin() << "Solving for" << this->f_.get_name();
        decltype(auto) expr = boost::yap::as_expr(e);
        solve_impl(f_, id_, expr);
        clamp_in_range();
        this->log().info_end();
    }

    // get the current solution
    FieldType &get() { return f_; }

    // force an execution of the equation
    // use always previous timestep value
    // for ddt
    FieldType &estimate() {
        // dont advance in time for estimate
        this->estimate_ = true;
        this->execute();
        this->estimate_ = false;
        return f_;
    }

    // if limits are specified the equation will enforce bounds after calling
    // solve
    void init_limits() {
        std::string min_name = f_.get_name() + "_min";
        std::string max_name = f_.get_name() + "_max";

        if (parameter_[min_name] || parameter_[max_name]) bounded_ = true;

        bounds_.first =
            read_or_default_coeff(min_name, std::numeric_limits<Scalar>::max());
        bounds_.second =
            read_or_default_coeff(max_name, std::numeric_limits<Scalar>::min());
    }

    // restricts field values to be in given range
    void clamp_in_range() {
        if (bounded_) {
            clamp_field_in_range(f_, bounds_.first, bounds_.second);
        }
    }
};

template <class FieldGradientType, template <typename> class SumABdWType>
class FieldGradientEquation : public FieldEquationBase<FieldGradientType> {

    // TODO use two separate fields for gradients of the Kernel
    // in the symmetric case the neighbour owner gradient field is just a
    // reference other wise the two fields are different fields

  protected:
    KernelGradientField &dWdx_;
    KernelGradientField &dWdxn_;

    SumABdWType<FieldGradientType> sum_AB_dW_s;

  public:
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
          dWdx_(objReg.get_object<KernelGradientField>("KerneldWdx")),
          dWdxn_(objReg.get_object<KernelGradientField>("KerneldWdxNeighbour")),
          sum_AB_dW_s(SumABdWType<FieldGradientType>(
              this->f_, this->np_, dWdx_, dWdxn_)) {}

    template <class RHS> void solve(RHS rhs, bool reset = true) {
        this->log().info_begin() << "Solving for " << this->f_.get_name();
        if (reset)
            solve_impl_res(this->f_, this->id_, rhs);
        else
            solve_impl(this->f_, this->id_, rhs);

        sum_AB_dW_s.a = 0;
        sum_AB_dW_s.ab = 0;
        this->clamp_in_range();
        this->log().info_end();
    }
};

template <
    class FieldType,
    template <typename>
    class SumType,
    template <typename>
    class SumABdWType>
class FieldValueEquation : public FieldEquationBase<FieldType> {

  protected:
    ScalarField &W_;

    KernelGradientField &dWdx_;

    KernelGradientField &dWdxn_;

    SumABdWType<FieldType> sum_AB_dW_s;

  public:
    typedef FieldType value_type;

    FieldValueEquation(
        const std::string &eqn_base_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        FieldType &f,
        bool symmetric = true)
        : FieldEquationBase<FieldType>(eqn_base_name, parameter, objReg, f),
          W_(objReg.get_object<ScalarField>("KernelW")),
          dWdx_(objReg.get_object<KernelGradientField>("KerneldWdx")),
          dWdxn_(objReg.get_object<KernelGradientField>("KerneldWdxNeighbour")),
          sum_AB_dW_s(
              SumABdWType<FieldType>(this->f_, this->np_, dWdx_, dWdxn_)) {};

    Sum_AB_sym<FieldType> sum_AB() {
        return Sum_AB_sym<FieldType>(this->f_, this->np_, this->id_);
    }

    Sum_AB_asym<FieldType> sum_AB_a() {
        return Sum_AB_asym<FieldType>(this->f_, this->np_);
    }

    Sum_AB_dW_asym<FieldType> sum_AB_dW_asym() {
        return Sum_AB_dW_asym<FieldType>(this->f_, this->np_, dWdx_, dWdxn_);
    }

    template <class RHS> void solve(RHS rhs, bool reset = true) {
        this->log().info_begin() << "Solving for " << this->f_.get_name();
        if (reset)
            solve_impl_res(this->f_, this->id_, rhs);
        else
            solve_impl(this->f_, this->id_, rhs);
        sum_AB_dW_s.a = 0;
        sum_AB_dW_s.ab = 0;
        this->clamp_in_range();
        this->log().info_end();
    }
};

using ScalarFieldEquation =
    FieldValueEquation<ScalarField, Sum_AB_sym, Sum_AB_dW_sym>;
using VectorFieldEquation =
    FieldValueEquation<VectorField, Sum_AB_sym, Sum_AB_dW_sym>;
using VectorFieldEquationA =
    FieldValueEquation<VectorField, Sum_AB_asym, Sum_AB_dW_asym>;

using ScalarGradientEquation =
    FieldGradientEquation<VectorField, Sum_AB_dW_sym>;
// TODO clarify why this is not a tensor field
using VectorGradientEquation =
    FieldGradientEquation<VectorField, Sum_AB_dW_sym>;

// using ScalarIntegralEquation = FieldGradientEquation<ScalarField,
// Sum_AB_dW_sym>;
#endif
