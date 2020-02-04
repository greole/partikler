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
        Scalar ret;

        Scalar h2 = h * h;
        Scalar h3 = h2 * h;
        Scalar h6 = h3 * h3;
        Scalar h9 = h6 * h3;

        Scalar pi = (Scalar)M_PI;

        if (2.0 * r > h && r <= h) {
            ret = (h - r) * (h - r) * (h - r) * r * r * r;

            return 32.0 / pi * h9 * ret;
        }
        if (r > 0 && 2 * r <= h) {
            ret = 2.0 * (h - r) * (h - r) * (h - r) * r * r * r - h6 / 64.0;

            return 32.0 / pi * h9 * ret;
        }

        return 0.0;
    };

    Scalar h;
};

CohesionAkinci::CohesionAkinci(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "fCohesion", {}, {"fCx", "fCy", "fCz"})),
      conti_(objReg.get_or_create_model<Conti>("Conti", parameter, objReg)),
      mp_(objReg.get_object<Generic<Scalar>>("specific_particle_mass")()),
      gamma_(read_or_default_coeff<Scalar>("gamma", 1.0)),
      h_(read_or_default_coeff<Scalar>("h", 1.0)),
      rho_0_(read_or_default_coeff<Scalar>("rho", 1.0)),
      n_(objReg.create_field<VectorField>("normal", {}, {"nx", "ny", "nz"})),
      pos_(objReg.get_object<VectorField>("Pos")) {}

void CohesionAkinci::execute() {

    log().info_begin() << "Computing surface normals";
    ScalarField &rho = conti_.get(time_.get_current_timestep());

    auto &dW = get_objReg().get_object<KernelGradientField>("KerneldWdx");
    sum_AB_dW_res_impl(n_, np_, dW, mp_ / rho.b());
    log().info_end();

    log().info_begin() << "Computing coehesive forces";
    iteration_ = time_.get_current_timestep();

    Scalar gamma = -1.0 * gamma_;

    auto norm = boost::yap::make_terminal(Norm_Wrapper());
    auto C = boost::yap::make_terminal(SplineFunc(h_));

    sum_AB(
        2.0 * rho_0_/(rho.a()+rho.b())*(
        mp_ * mp_ * gamma * C(norm(ab(pos_))) * ab(pos_) / norm(ab(pos_)) +
        mp_ * gamma * ab(n_))
            );

    log().info_end();

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(SURFACETENSION, Akinci);
