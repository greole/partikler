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

#ifndef CGALHelper_H
#define CGALHelper_H

#include "CGALTYPEDEFS.hpp" // for Facet_handle, CGALVector
#include "Helper.hpp"
#include "Vec3.hpp" // for Vec3

// A collection of CGAL related functions

CGALVector normalise(const CGALVector &v);

CGALVector vectorRejection(const CGALVector &a, const CGALVector &b);

bool approachesEdge(
    Point A,      // Beginning edge
    Point B,      // End edge
    Point P,      // Particle start position
    CGALVector &S // Slide vector
);

CGALVector facet_normal(const Facet &f);

CGALVector surfaceVector(const Facet &f);

struct Point3D {
    // 3*4bytes = 12bytes
    float x, y, z;
};

struct HitPoint {
    bool hit;
    Point X;
};

struct STLSurfaceDist {
    // Stores the distance of particles on different
    // STL surfaces

    float len;
    CGALVector on;
    CGALVector no;
};

HitPoint approximateEdgeHit(
    Point &A,      // Beginning edge
    Point &B,      // End edge
    Point &P,      // Particle start position
    CGALVector &S, // Slide vector
    float tolerance);

struct Matrix {
    // x0 y0 z0
    // x1 y1 z1
    // x2 y2 z2
    CGALVector x;
    CGALVector y;
    CGALVector z;
};

CGALVector surfaceProject(CGALVector N, CGALVector S);

Matrix mult(Matrix a, Matrix b);

Matrix inv(Matrix X);

CGALVector rotate(Matrix a, CGALVector b);

Matrix facetCoordSys(Facet f);

Matrix createRotationMatrix(Facet a, Facet b);

double length(CGALVector a);

Matrix createRotationMatrix(CGALVector N1, CGALVector N2);

struct CommonVertices {
    int n;
    Point A;
    Point B;
};

CommonVertices commonVertices(const Facet &F1, const Facet &F2);

// TODO compute face2facePath of type std::vector<Facet> from
// (Point P1, Point P2, Facet F1, Facet F2) {
// Step 1. check if F1 and F2 have two common vertices
// if so stop
// otherwise compute face2facePath of neighbour of F1
// if already recursed then dont do another recursion
// if no two common points try next facet neighbour
typedef std::vector<Facet_handle> Path;

std::pair<bool, Path> searchPath(Path p);

Point projectedPoint(Path p, Point P);

Vec3 rotate(Matrix R, Vec3 a);

Triangle facetToTriangle(const Facet &facet);

// Functors
struct Compute_Facet_Area {
    double operator()(const Facet &f) const;
};

struct Compute_Facet_Normal {
    inline CGALVector operator()(const Facet &f) const {
        return facet_normal(f);
    }
};

struct Compute_Facet_Directions {
    inline CGALVector operator()(const Facet &f) const {
        HalfedgeConstHandle h = f.halfedge();
        Point p1 = h->vertex()->point();
        Point p2 = h->next()->vertex()->point();
        CGALVector n = p2 - p1;
        return n / std::sqrt(n * n);
    }
};

struct movePoint {
    CGALVector n;
    float x, y, z;

    inline movePoint(CGALVector normal, float dx) {
        n = normal;
        x = dx;
        y = dx;
        z = dx;
    }

    inline Point operator()(const Point &p) const {
        return Point(p.x() + x * n.x(), p.y() + y * n.y(), p.z() + z * n.z());
    }
};

struct EdgeNormal {
    CGALVector EN;
    CGALVector FN;
};

EdgeNormal inwardPointingEdgeNormal(Point &A, Point &B, Facet &f);

Point limited_project_point(
    Point P,
    float dx_,
    float normal_dist,
    float angle,
    CGALVector CP,
    double max_lambda);

struct Generate_Points_at_Facets {

    const float dx_;
    // const float dx2_;
    std::vector<Point> &points_;
    // std::vector<size_t>& number_points_facet_;
    std::vector<Facet_handle> &initial_facets_;
    // std::vector<CGALVector>& facet_vectors_;

    size_t facetId;

    Generate_Points_at_Facets(
        const float dx,
        std::vector<Point> &points,
        std::vector<Facet_handle> &initial_facets)
        : dx_(dx),
          // dx2_(dx*dx),
          points_(points),
          // number_points_facet_(number_points_facet),
          initial_facets_(initial_facets),
          // facet_vectors_(facet_vectors),
          facetId(0) {}

    void operator()(Facet &f);
};

std::ostream &operator<<(std::ostream &os, Point const &p);

STLSurfaceDist compute_STLSurface_dist(
    Point opos, Point npos, Facet_handle start, Facet_handle end);

float operator*(Vec3 &x, CGALVector &y);

#endif
