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

#include "Solenthaler.hpp"

#include "Time.hpp"


SolenthalerGradient::SolenthalerGradient(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : ScalarGradientEquation(
        "PressureGradient",
        parameter,
        objReg,
        objReg.get_object<ScalarField>("p")),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      pressure_(objReg.get_object<ScalarFieldEquation>("Pressure")),
      vel_(objReg.get_object<VectorFieldEquation>("Momentum")),
      pos_(objReg.get_object<VectorFieldEquation>("Position")),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      eta_(read_or_default_coeff<Scalar>("eta", 1.0)),
      mp_(objReg.get_object<Generic<Scalar>>("specific_particle_mass")())
{}


void SolenthalerGradient::execute() {


    // Check convergence criteria

    for(int iter=0; iter<3;iter++) {

        // Predict velocity
        int it = time_.get_current_timestep();

        auto & vel = vel_.get(it + iter);

        auto & pos = pos_.get(it + iter);

        auto & rho = conti_.get(it + iter);

        // TODO compute density variation

        auto & p = pressure_.get(it + iter);

        // update pressure

        Scalar beta = 2.0 * time_.get_deltaT() * mp_ * mp_ / rho_0_ / rho_0_;

        auto &dW = this->get_objReg().template get_object<KernelGradientField>(
            "KerneldWdx");


        // TODO implement a lazy variant
        VectorField sumABone(f_.size(), {0.0, 0.0, 0.0});
        ScalarField sumABtwo(f_.size(), 0.0);
        for(size_t ab=0; ab<np_.size();ab++) {
            size_t i = np_[ab].ownId;
            size_t j = np_[ab].neighId;

            Scalar sumABtwores = dW[ab]*dW[ab];

            sumABone[i] += dW[ab];
            sumABone[j] -= dW[ab];

            sumABtwo[i] += sumABtwores;
            sumABtwo[j] -= sumABtwores;

        }

        ScalarField ptilda(f_.size(), 0.0);
        Scalar rho_err = 0.0;

        solve_impl(
            ptilda,
            -(rho - rho_0_) / (beta * (-(sumABone * sumABone) - sumABtwo)));

        for(size_t i=0; i<f_.size(); i++) {
            p[i] += ptilda[i];
        }

        auto sum_AB_e = Sum_AB_dW_sym<VectorField>(f_, np_, dW, dW);
        auto sum_AB_dW_e = boost::yap::make_terminal(sum_AB_e);

        solve(
            mp_*mp_*sum_AB_dW_e(
                p.a()/(rho.a()*rho.a())
              + p.b()/(rho.b()*rho.b())
                ),
            true
            );
    }

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(PRESSURE, SolenthalerGradient);
