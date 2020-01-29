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

#include "Models.hpp"
#include "Field.hpp"            // for ScalarField, KernelGradientField
#include "FieldOps.hpp"         // for sum_AB_impl
#include "Time.hpp"
#include "Logger.hpp"           // for MSG, Logger
#include "Object.hpp"           // for SPHObject, ModelType
#include "ObjectRegistry.hpp"   // for ObjectRegistry
#include "SearchCubes.hpp"      // for NeighbourFieldAB
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

template <class T, class U> class FieldEquation : public Model {

  protected:
    TimeGraph &time_;

    NeighbourFieldAB &np_;

    ScalarField &W_;

    KernelGradientField &dW_;

    T &f_;

    // decltype(gradient_type<T>::type) &ddxf_; // TODO compute type
    U &df_;

    // store previous results
    std::map<int, T> prev_;

    // the current iteraton
    int iteration_;

  public:
    FieldEquation(
        const std::string &field_name,
        YAML::Node parameter,
        ObjectRegistry &objReg,
        T &f)
        : Model(field_name, parameter, objReg),
          time_(objReg.get_object<TimeGraph>("TimeGraph")),
          np_(objReg.get_object<NeighbourFieldAB>("neighbour_pairs")),
          W_(objReg.get_object<ScalarField>("KernelW")),
          dW_(objReg.get_object<KernelGradientField>("KerneldWdx")), f_(f),
          df_(objReg.create_field<VectorField>(
              "d" + f_.get_name(),
              zero<VectorField::value_type>::val,
              {"dx" + f_.get_name(),
               "dy" + f_.get_name(),
               "dz" + f_.get_name()})) {};

    // get result for iteration i
    // if result is not cached solve gets executed

    ScalarField &W() { return W_; };

    NeighbourFieldAB &N() { return np_; };

    KernelGradientField &dWdx() { return dW_; };

    ScalarField &sum_AB(float particle_mass) {
        sum_AB_impl(particle_mass, f_, np_, W_);
        return f_;
    };

    template <class RHS> ScalarField &sum_AB(RHS rhs) {
        sum_AB_impl(f_, np_, rhs * W_);
        return f_;
    }

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
    T &get(int iteration) {
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

    U &get_dx(int iteration) {
        if (iteration == iteration_) {
            std::cout << get_name()
                      << "DEBUG using cached solution for iteration "
                      << iteration << std::endl;
            return df_;
        } else {
            this->execute();
            return df_;
        }
    }
};

using ScalarFieldEquation = FieldEquation<ScalarField, VectorField>;
using VectorFieldEquation = FieldEquation<VectorField, VectorField>;

#endif
