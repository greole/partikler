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

#ifndef PARTIKLER_VEC3_INCLUDED
#define PARTIKLER_VEC3_INCLUDED

#include "Scalar.hpp"
#include <array>
#include <iostream>
#include <vector>

struct Vec3 : std::array<Scalar, 3> {

    Vec3 &operator=(const Vec3 &x);

    Vec3 &operator+=(const Vec3 &x);

    Vec3 &operator-=(const Vec3 &x);

    Vec3 &operator*=(const float x);
};

// Eager Operators
//
// Here are the operators on the Vec3 type
// which are eager, since they are only used as
// inner type for the lazy field functions

// scalar multiplication
Vec3 operator*(Scalar const a, Vec3 const &x);
Vec3 operator*(Vec3 const &x, Scalar const a);

// scalar division
Vec3 operator/(const Vec3 &x, Scalar a);

// dot product
float operator*(const Vec3 &x, const Vec3 &y);

// addition
Vec3 operator+(const Vec3 &x, const Vec3 &y);

Vec3 operator-(const Vec3 &x, const Vec3 &y);

std::ostream &operator<<(std::ostream &os, Vec3 const &f);

// A pair of Vec3 to store particle particle  distances on stl surfaces
struct VectorPair {
    Vec3 on;
    Vec3 no;
};

Scalar squared_length(Vec3 v);

void translatePoints(std::vector<Vec3> &points, Vec3 translate);

void scalePoints(std::vector<Vec3> &points, Vec3 scale);

void scalePoints(std::vector<Vec3> &points, Scalar scale);

Scalar mag(Vec3 v);

#endif
