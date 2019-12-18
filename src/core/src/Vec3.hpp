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

#include <array>
#include <iostream>

struct Vec3: std::array<float, 3> {

    Vec3& operator=(const Vec3 &x);

};

// Eager Operators
//
// Here are the operators on the Vec3 type
// which are eager, since they are only used as
// inner type for the lazy field functions

// scalar multiplication
Vec3 operator*(float a, const Vec3& x);

// scalar division
Vec3 operator/(const Vec3& x, float a);

// dot product
float operator*(const Vec3& x, const Vec3& y);

// addition
Vec3 operator+(const Vec3& x, const Vec3& y);

Vec3 operator-(const Vec3& x, const Vec3& y);

std::ostream &operator<<(std::ostream &os, Vec3 const &f);

// A pair of Vec3 to store particle particle  distances on stl surfaces
struct VectorPair {
    Vec3 on;
    Vec3 no;
};

#endif
