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

#include "Akinci.hpp"

#include "Conti.hpp"
#include "Time.hpp"
#include <math.h>

struct SplineFunc {
    SplineFunc(Scalar i) { h = i; };
    Scalar operator()(Scalar r) {
        // TODO refactor
        Scalar h2 = h * h;
        Scalar h3 = h2 * h;
        Scalar h6 = h3 * h3;
        Scalar h9 = h6 * h3;

        Scalar r3 = r * r * r;
        Scalar hr3 = (h - r) * (h - r) * (h - r);

        Scalar pi = (Scalar)M_PI;

        Scalar c = 32.0 / pi / h9;
        Scalar hr3r3 = hr3 * r3;

        if (2.0 * r > h && r <= h) {
            return c * hr3r3;
        }
        if (r > 0 && 2 * r <= h) {
            return c * (2.0 * hr3r3 - h6 / 64.0);
        }

        return 0.0;
    };

    Scalar h;
};

Akinci::Akinci(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          "Cohesion",
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "f", {}, {"fCoex", "fCoey", "fCoez"})),
      conti_(objReg.get_object<ScalarFieldEquation>("Conti")),
      mp_(objReg.get_object<ScalarField>("mp")),
      gamma_(read_or_default_coeff<Scalar>("gamma", 1.0)),
      h_(read_or_default_coeff<Scalar>("h", 1.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho_0", 1.0)),
      pos_(objReg.get_object<VectorField>("Pos")) {}

void Akinci::execute() {

    log().info_begin() << "Computing coehesive forces";
    iteration_ = time_.get_current_timestep();

    Scalar gamma = -1.0 * gamma_;

    ScalarField &rho = conti_.get();

    auto saba = sum_AB_a();
    auto sum_AB_e = boost::yap::make_terminal(saba);
    auto norm_e = Norm_Wrapper();
    auto norm = boost::yap::make_terminal(norm_e);
    auto spline_f = SplineFunc(h_);
    auto C = boost::yap::make_terminal(spline_f);

    // TODO do this via a solve method
    VectorFieldAB dist(this->np_.size(), {0, 0, 0});
    solve_inner_impl(this->np_, dist, ab(pos_));

    Scalar c = 2.0 * rho_0_ *  gamma;

    solve(sum_AB_e(
              mp_.a() * mp_.b() / (rho.a() + rho.b()) * (C(norm(dist)) * ab(pos_) / norm(dist))));

    log().info_end();

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(SURFACETENSION, Akinci);
