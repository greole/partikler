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

#include "CGALHelper.hpp"

CGALVector normalise(const CGALVector& v) {
    // normalise the vector v
    return v/std::sqrt(v.squared_length());
}

CGALVector vectorRejection(const CGALVector & a, const CGALVector & b) {
    return normalise(b - (a*b)*a);
}

bool approachesEdge(
    Point A,      // Beginning edge
    Point B,      // End edge
    Point P,      // Particle start position
    CGALVector &S // Slide vector
    ) {
    // Check if Vector S from point P approaches an Edge
    // Does not check if it hits the edge
    Line AB(A, B);

    // Point N = P + S;

    // Point XN is the normal intersection Point
    Point X = AB.projection(P);

    // distance from intersection point X to N
    CGALVector PX = X - P;

    // Check if slide vector and edge normal
    // have same direction, ie particle move towards
    // edge
    return PX * S > 0;
}

CGALVector facet_normal(const Facet& f) {
    // Compute the normal vector of the facet f

    HalfedgeConstHandle h = f.halfedge();
    Point p1 = h->vertex()->point();
    Point p2 = h->next()->vertex()->point();
    Point p3 = h->next()->next()->vertex()->point();

    CGALVector n = CGAL::cross_product(p2 - p1, p3 - p1);

    // std::cout
    //   << "[DEBUG SLIDE]"
    //   << " P1 " << p1
    //   << " P2 " << p2
    //   << " P3 " << p3
    //   << std::endl;

    return n / std::sqrt(n * n);
}


CGALVector surfaceVector(const Facet& f) {
    // Compute a vector on the facet f

    HalfedgeConstHandle h = f.halfedge();
    Point p1 = h->vertex()->point();
    Point p2 = h->next()->vertex()->point();

    return normalise(p2-p1);
}

HitPoint approximateEdgeHit(
    Point &A,      // Beginning edge
    Point &B,      // End edge
    Point &P,      // Particle start position
    CGALVector &S, // Slide vector
    float tolerance
) {
    // Given a tolerance this function checks if the Vector S from point P hits
    // the edge AB
  Line G1 = Line(A, B);
  Line G2 = Line(P, P + S);
  auto res = intersection(G1, G2);
  Point X;
  bool intersects = assign(X, res);
  if (intersects) return {intersects, X};

  // if intersects was false, it might be due to floating point precision
  // assume a hit since the slide vector is a surface vector
  CGALVector D1 = B - A;
  CGALVector D2 = S;
  CGALVector N2 = cross_product(D2, cross_product(D1, D2));
  float t = ((P - A) * N2) / (D1 * N2);
  X = A + t * D1;
  if ((X-P).squared_length() < S.squared_length()) return {true, X};
  return {false, X};



  // // TODO REVISE
  // // At this point it is already known that starting from P S approaches AB
  // // Hence test if S would hit AB within lambda 0..1 and
  // // and if length |PX| < |S| then X is A + lambda AB

  //   Point OO {0,0,0};
  //   CGALVector BO = B - OO;
  //   CGALVector PS = (P + S) - OO;
  //   CGALVector BB = CGAL::cross_product(BO, PS);
  //   float distance = (A - P) * BB / std::sqrt(BB.squared_length());

  //   if (distance * distance < tolerance) {
  //       intersects = true;
  //       CGALVector D1 = B - A;
  //       CGALVector D2 = S;
  //       CGALVector N2 = cross_product(D2, cross_product(D1, D2));
  //       float t = ((P - A) * N2) / (D1 * N2);
  //       X = A + t * D1;

  //       // std::cout
  //       //   << "[DEBUG] approximateEdgeHit "
  //       //   << " distance " << distance
  //       //   << " t " << t
  //       //   << " A " << A
  //       //   << " B " << B
  //       //   << " P " << B
  //       //   << " D1 " << D1
  //       //   << " S " << S
  //       //   << " X " << X
  //       //   << std::endl;
  //   }
  //   return {intersects, X};
}

CGALVector surfaceProject(CGALVector N, CGALVector S) {
    return S - (N * S) * N;
}

Matrix mult(Matrix a, Matrix b) {

    return {CGALVector(
            a.x[0] * b.x[0] + a.y[0] * b.x[1] + a.z[0] * b.x[2],
            a.x[1] * b.x[0] + a.y[1] * b.x[1] + a.z[1] * b.x[2],
            a.x[2] * b.x[0] + a.y[2] * b.x[1] + a.z[2] * b.x[2]),
            CGALVector(
                a.x[0] * b.y[0] + a.y[0] * b.y[1] + a.z[0] * b.y[2],
                a.x[1] * b.y[0] + a.y[1] * b.y[1] + a.z[1] * b.y[2],
                a.x[2] * b.y[0] + a.y[2] * b.y[1] + a.z[2] * b.y[2]),
            CGALVector(
                a.x[0] * b.z[0] + a.y[0] * b.z[1] + a.z[0] * b.z[2],
                a.x[1] * b.z[0] + a.y[1] * b.z[1] + a.z[1] * b.z[2],
                a.x[2] * b.z[0] + a.y[2] * b.z[1] + a.z[2] * b.z[2])};
}



Matrix inv(Matrix X) {
    // a b c
    // d e f
    // g h i
    double a = X.x[0];
    double b = X.y[0];
    double c = X.z[0];

    double d = X.x[1];
    double e = X.y[1];
    double f = X.z[1];

    double g = X.x[2];
    double h = X.y[2];
    double i = X.z[2];

    double A =   e*i-f*h;
    double B = -(d*i - f*g);
    double C =   d*h-e*g;

    double D = -(b*i-c*h);
    double E =   a*i-c*g;
    double F = -(a*h-b*g);

    double G =   (b*f-c*e);
    double H =  -(a*f-c*d);
    double I =   (a*e-b*d);

    // Rule of Sarrus
    double id = 1./(a*A+b*B+c*C);

    return {id*CGALVector(A, B, C),
            id*CGALVector(D, E, F),
            id*CGALVector(G, H, I)};
}

EdgeNormal inwardPointingEdgeNormal(
                                    Point &A,
                                    Point &B,
                                    Facet &f
                                    ) {
  // compute the edge normal which points inwards
  Line AB {A, B};
  CGALVector ON =  facet_normal(f);

  // inwards pointing vector
  // assume that other component points always outwards

  // TODO refactor
  // Find the other, opposing vertex
  auto e = f.halfedge();
  Point nC = e->opposite()->vertex()->point();
  if (nC == A || nC == B) {
    nC = e->next()->vertex()->point();
  }
  if (nC == A || nC == B) {
    nC = e->vertex()->point();
  }

  CGALVector NAB = CGAL::cross_product(B-A, ON);
  NAB = NAB/std::sqrt(NAB.squared_length());
  Point ABe = AB.point(0);
  CGALVector NC = nC-ABe;
  // std::cout
  //   << "[DEBUG SLIDE] inwardPointingEdgeNormal"
  //   << " NAB" << NAB
  //   << " AB " << AB
  //   << " nC " << nC
  //   << " ON " << ON
  //   << " AB.point(0) " << ABe
  //   << std::endl;

  if (NAB * NC < 0) {
    NAB = -NAB;
  }

  return {NAB, ON};
}

CGALVector rotate(Matrix a, CGALVector b)
{
    return {
        a.x[0] * b.x() + a.y[0] * b.y() + a.z[0] * b.z(),
        a.x[1] * b.x() + a.y[1] * b.y() + a.z[1] * b.z(),
        a.x[2] * b.x() + a.y[2] * b.y() + a.z[2] * b.z()
    };
}

Matrix facetCoordSys(Facet f) {
    CGALVector norm = facet_normal(f);
    Point p1 = f.facet_begin()->vertex()->point();
    Point p2 = f.facet_begin()->next()->vertex()->point();
    CGALVector surf = p2 - p1;
    surf = surf/std::sqrt(surf.squared_length());
    CGALVector third = cross_product(surf, norm);

    return { surf, norm, third };
}

Matrix createRotationMatrix(Facet a, Facet b) {
    // create a rotation matrix from two rotated facets
    // uses e2*e1^-1 as rotation matrix
    Matrix ea = facetCoordSys(a);    // unit vectors of new facet
    Matrix eb = facetCoordSys(b); // unit vectors of old facet
    Matrix ieb = inv(eb);

    // std::cout
    //   << "[DEBUG SLIDE] create rotation matrix\n"
    //   << "[DEBUG SLIDE]"
    //   << " ea.x " << ea.x
    //   << " ea.y " << ea.y
    //   << " ea.z " << ea.z
    //   << " eb.x " << eb.x
    //   << " eb.y " << eb.y
    //   << " eb.z " << eb.z
    //   << " ieb.x " << ieb.x
    //   << " ieb.y " << ieb.y
    //   << " ieb.z " << ieb.z
    //   << std::endl;

    return mult(ea, ieb);
}

double length(CGALVector a) {
    return std::sqrt(a.squared_length());
}

Matrix createRotationMatrix(CGALVector N1, CGALVector N2) {
    //  https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d  
    double x0 = N1*N2;

    // handle no rotation case
    if (x0 == 1) {
        return {
            CGALVector(1, 0, 0),
            CGALVector(0, 1, 0),
            CGALVector(0, 0, 1)};

    }
    double x1 = length(cross_product(N1,N2));
    double x2 = 0;

    Matrix G = {
        CGALVector(x0, x1, x2),
        CGALVector(-x1, x0, x2),
        CGALVector(0, 0, 1)};

    Matrix F = { N1, vectorRejection(N1, N2), cross_product(N1, N2) };
    Matrix invF = inv(F);

    return mult(F, mult(G, invF));
}

CommonVertices commonVertices(const Facet &F1, const Facet &F2)
{
    // find a pair of common vertices for two adjacent facets

    Point P; // Pivot Point
    Point T; // Test Point
    Point retA;

    bool found = false;
    int n = 0;

    auto f1_begin = Facet(F1).facet_begin();
    auto f1 = F1.facet_begin();

    auto f2_begin = Facet(F2).facet_begin();

    do {
        // Set pivot point;
        P = f1->vertex()->point();
        auto f2 = F2.facet_begin();
        do {
            T = f2->vertex()->point();
            // std::cout << "[FACE2FACEDIST]"
            //           << " P " << P
            //           << " T " << T
            //           << std::endl;
            if (T == P) {
              n++;

              // std::cout << "[FACE2FACEDIST]"
              //           << " FOUND " << std::endl;
                // already found a point pair
                if (found) return {n, retA, T};

                // found first point
                found = true;
                retA = T;
            }
            f2++;
        } while (f2 != f2_begin);
        f1++;
    } while (f1 != f1_begin);

    return {n, retA, T};
}

Point limited_project_point(
    Point P,
    float dx_,
    float normal_dist,
    float angle,
    CGALVector CP,
    double max_lambda) {
    // Check if Deltax yields a lambda > max_lambda
    // TODO refactor so that this check is performed only once
    const float initial_lambda = dx_ / angle * 0.5;
    float lambda;
    if (initial_lambda > max_lambda) {
        lambda = max_lambda;
    } else {
        lambda = normal_dist / angle;
    }

    return P + lambda * CP;
}

void Generate_Points_at_Facets::operator()(Facet &f) {

    Point A = f.halfedge()->vertex()->point();
    Point B = f.halfedge()->next()->vertex()->point();
    Point C = f.halfedge()->opposite()->vertex()->point();

    // const CGALVector facet_normal = Compute_Facet_Normal()(f);
    const CGALVector edge_normal = inwardPointingEdgeNormal(A, B, f).EN;
    const CGALVector CA = C - A;
    const CGALVector CB = C - B;

    float dist_base = dx_ / 2.0; // Cartesian distance from base edge
    float dist_points =
        dx_ / 2.0; // Cartesian dist from beginning of point line

    // Unit vectors of triangle sides
    CGALVector UCA = normalise(CA);
    CGALVector UCB = normalise(CB);
    CGALVector UCE = normalise(edge_normal);

    // length of CA divided by length of UCA divided by two
    float max_lambda_s = length(CA) / 2.0;
    float max_lambda_e = length(CB) / 2.0;

    float angle_s = UCA * UCE;
    float angle_e = UCB * UCE;

    float dx2 = dx_*dx_;

    int line_ctr = 0;
    srand(time(NULL));
    float first_dist = dx2/4.0;

    do {
        Point S = limited_project_point(
            A, dx_, dist_base, angle_s, UCA, max_lambda_s);
        Point E = limited_project_point(
            B, dx_, dist_base, angle_e, UCB, max_lambda_e);

        const CGALVector SE = E - S;
        const float lenSEsqr =  SE.squared_length();

        const CGALVector USE = normalise(SE);

            do {

                // TODO if initial point is already outside of SE throw a dice
                // an set it on SE/2
                if ( first_dist > lenSEsqr) {
                    float ratio = lenSEsqr/first_dist;

                    float  dice = ((float)(rand() % 100))/100.0;
                    if ( dice < ratio ) {

                    Point P = S + 0.5*SE;

                    points_.push_back(P);
                    initial_facets_.push_back(&f);

                    }

                break;
            }

            Point P = S + USE * dist_points;
            dist_points += 1.5 * dx_;

            points_.push_back(P);
            initial_facets_.push_back(&f);
        } while (
            // dist_points is used to avoid sqrt(|SE|)
            dist_points * dist_points < lenSEsqr - dx2
            );

        if (line_ctr % 2 == 0)
            dist_points = dx_ / 2.0;
        else
            // Even start points along facet edge need to be pushed slightly
            // inwards to avoid edge detection misses
            dist_points = dx_ * 0.001;

        line_ctr++;
        dist_base += dx_ / 2.0;

    } while (dist_base * dist_base <
             (angle_s * CA).squared_length() - (dx_ * dx_));

    facetId++;
}

std::pair<bool, Path> searchPath(Path p)
{
  // return empty path if search failed
  int nhops = p.size();

  Facet f1 = Facet(*p[nhops-2]);
  Facet f2 = Facet(*p[nhops-1]);

  CommonVertices commonVerts = commonVertices(f1, f2);
  // success if last facets are connected
  if (commonVerts.n == 2) return {true, p};

  if (commonVerts.n == 1) {
    // pick a neighbour 

    auto h = f1.facet_begin();

    do {
      // Check if opposite halfedge exists
      if (h == 0) {
          h++;
          continue;
      }

      // get an edge an get opposite facet
      Facet_handle next_facet_handle = h->opposite()->facet();
      Facet next_facet = Facet(*next_facet_handle);

      // check if current neighbour facet has two common vertices
      commonVerts = commonVertices(next_facet, f2);
      if (commonVerts.n == 2) {
        p.insert(p.begin()+nhops-1, next_facet_handle);
        return {true, p};
      }
      h++;

    } while ( h != f1.facet_begin());
  }

  return {false, p};
}

Point projectedPoint(Path p, Point P)
{
    int nhops = p.size();
    // if path is empty just return P. Results in Cartesian distance
    if (nhops < 2) return P;

    Facet f1 = Facet(*p[nhops-2]);
    Facet f2 = Facet(*p[nhops-1]);

    CGALVector N1 = facet_normal(f1);
    CGALVector N2 = facet_normal(f2);

    if (N1 == N2) {
        p.pop_back();
        return projectedPoint(p, P);
    }

    Matrix RO = createRotationMatrix(N2, N1);

    // RO has its origin on the connecting edge, thus a transformation of
    // the coordinate system is needed
    CommonVertices PAB = commonVertices(f1, f2);
    Point A = PAB.A;
    Point B = PAB.B;
    CGALVector AB = B - A;
    CGALVector AP = P - A;

    CGALVector ABP1 = rotate(RO, AP);

    p.pop_back();
    return projectedPoint(p, A + ABP1);
}


Vec3 rotate(Matrix R, Vec3 a) {
  CGALVector v {a[0], a[1], a[2]};
  CGALVector ret {R.x*v, R.y*v, R.z*v};
  return {(float)ret[0], (float)ret[1], (float)ret[2]};
}

Triangle facetToTriangle(const Facet& facet) {
  // Needed for the random point generator
    return Triangle(
           facet.halfedge()->vertex()->point(),
           facet.halfedge()->next()->vertex()->point(),
           facet.halfedge()->opposite()->vertex()->point()
       );
}


double Compute_Facet_Area::operator()(const Facet& f) const {
        return K::Compute_area_3()(
            f.halfedge()->vertex()->point(),
            f.halfedge()->next()->vertex()->point(),
            f.halfedge()->opposite()->vertex()->point() );
}


STLSurfaceDist compute_STLSurface_dist(
    Point opos, Point npos,
    Facet_handle start, Facet_handle end) {

    if (start == end) {
        // if particles on same facet return Cartesian distance
        CGALVector lenVo = npos - opos;
        CGALVector lenVn = opos - npos;
        return {(float) length(lenVo), lenVo, lenVn};
    }

    std::pair<bool, Path> ret_path = searchPath({start, end});

    if (!ret_path.first) {
        // if path search failed use Cartesian distance 
        CGALVector lenVo = npos - opos;
        CGALVector lenVn = opos - npos;
        return {(float) length(lenVo), lenVo, lenVn};
    }

    Path path = ret_path.second;

    Point nposp, oposp;
    nposp = projectedPoint(path, npos);
    reverse(path.begin(), path.end());
    oposp = projectedPoint(path, opos);

    CGALVector lenVo = nposp - opos;
    CGALVector lenVn = oposp - npos;

    return {(float) length(lenVo), lenVo, lenVn};
}
