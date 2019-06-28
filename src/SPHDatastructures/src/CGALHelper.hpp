#ifndef CGALHelper_H
#define CGALHelper_H

#include "CGALTYPEDEFS.hpp"
#include "Helper.hpp"

CGALVector normalise(const CGALVector& v) {
  // normalise the vector v
  return v/std::sqrt(v.squared_length());
}

CGALVector vectorRejection(const CGALVector & a, const CGALVector & b) {
  return normalise(b - (a*b)*a);
};

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

struct HitPoint {
  bool hit;
  Point X;
};


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

struct Matrix {
  // x0 y0 z0
  // x1 y1 z1
  // x2 y2 z2
  CGALVector x;
  CGALVector y;
  CGALVector z;

};

CGALVector surfaceProject(CGALVector N, CGALVector S) {
  return S - (N * S) * N;
};

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

CGALVector rotate(Matrix a, CGALVector b) {
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

struct CommonVertices{
  int n;
  Point A;
  Point B;
};

CommonVertices commonVertices(const Facet &F1, const Facet &F2) {
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

// TODO compute face2facePath of type std::vector<Facet> from
// (Point P1, Point P2, Facet F1, Facet F2) {
// Step 1. check if F1 and F2 have two common vertices
// if so stop 
// otherwise compute face2facePath of neighbour of F1
// if already recursed then dont do another recursion
// if no two common points try next facet neighbour
typedef std::vector<Facet_handle> Path;

std::pair<bool, Path> searchPath(Path p) {
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


Point projectedPoint(Path p, Point P) {
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

// std::pair<CGALVector, CGALVector>
// // face2faceDistance(Point P1, Point P2, Facet F1, Facet F2) {
// face2faceDistance(Point P1, Point P2, Path path) {
//   // returns distance vector, ie P2-P1 but P2 has been rotated
//   // onto F1

//   CGALVector N1 = facet_normal(F1);
//   CGALVector N2 = facet_normal(F2);
//   CGALVector P1P2 = P2 - P1;
//   CGALVector P2P1 = P1 - P2;

//   if (N1 == N2)
//     {
//     // std::cout << "[FACE2FACEDIST]"
//     //           << " SKIP FACETS HAVE SAME NORMAL"
//     //           << std::endl;
//     return {P1P2, P2P1};
//     }

//   // A,B are the points which both facets have in common
//   CommonVertices PAB = commonVertices(F1, F2);
//   // std::cout << "[FACE2FACEDIST]"
//   //           << " FOUND " << PAB.n
//   //           << " COMMON VERTICES " << PAB.n
//             // << std::endl;

//   int nCommon = PAB.n;
//   CGALVector AB;
//   Point A = PAB.A;
//   Point B = PAB.B;
//   if (nCommon < 1) {
//     // std::cout << "[FACE2FACEDIST]"
//     //           << " returning Cartesian distance " 
//     //           << std::endl;
//     return {P1P2, P2P1};
//   }
//   if (nCommon == 2) {
//     A = PAB.A;
//     B = PAB.B;
//     AB = B-A;
//   }

//   if (nCommon == 1) {
//     return {P1P2, P2P1};
//     A = PAB.A;
//     AB = surfaceVector(F1);
//     B = A + AB;
//   }

//   // Rotate to owner face
//   Matrix RO = createRotationMatrix(N2, N1); //, AB);

//   // Rotate to neighbour face
//   Matrix RN = createRotationMatrix(N1, N2); // , -AB);

//   // Find intersection with edge
//   CGALVector P1P2P = surfaceProject(N1, P1P2);

//   // TODO invalid for one hit point
//   // A, B builds a line, which might be to short 
//   HitPoint x = approximateEdgeHit(A, B, P1, P1P2P, 1e-16);

//   // point to intersection connection vectors
//   CGALVector P2X = x.X - P2;
//   CGALVector P1X = x.X - P1;

//   // rotated connection vectors
//   CGALVector ABP1 = rotate(RN, P1-x.X);
//   CGALVector ABP2 = rotate(RO, P2-x.X);

//   // other side of the connection vector
//   Point P2P = x.X + ABP2;
//   Point P1P = x.X + ABP1;

//   CGALVector reto = P2P - P1;
//   CGALVector retn = P1P - P2;

//   // CGALVector foo = (rotate(RN, x.X-P1) - ( P2 - Point (0,0,0)));

//   // std::cout
//   //   << "[FACE2FACEDIST]"
//   //   << " N1 " << N1
//   //   << " N2 " << N2
//   //   << " P1P2P " << P1P2P
//   //   << " RN2 " << rotate(RO, N2)
//   //   << " RN1 " << rotate(RN, N1)
//   //   << " P1 " << P1
//   //   << " P2 " << P2
//   //   << " X " << x.X
//   //   << " P1X " << P1X
//   //   << " P2X " << P2X
//   //   << " ABP1 " << ABP1
//   //   << " ABP2 " << ABP2
//   //   << " P2P " << P2P
//   //   << " reto " << reto
//   //   << " RNretn1 " << -rotate(RO, reto)
//   //   << " RNretn2 " <<  rotate(RO, -reto)
//   //   << " RNretn3 " << -rotate(RN, reto)
//   //   << " RNretn3 " <<  rotate(RN, -reto)
//   //   << " RNretn4 " << foo
//   //   << " retn " << P1P - P2
//   //   << std::endl;

//   return {reto, retn};
// }


Vector rotate(Matrix R, Vector a) {
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

// Functors
struct Compute_Facet_Area
{
    double operator()(const Facet& f) const {
    return K::Compute_area_3()(
      f.halfedge()->vertex()->point(),
      f.halfedge()->next()->vertex()->point(),
      f.halfedge()->opposite()->vertex()->point() );
    }
};

struct Compute_Facet_Normal
{
    inline CGALVector operator() (const Facet& f) const {
      return facet_normal(f);
    }
};

struct Compute_Facet_Directions
{
    inline CGALVector operator() (const Facet& f) const {
        HalfedgeConstHandle h  = f.halfedge();
        Point              p1 = h->vertex()->point();
        Point              p2 = h->next()->vertex()->point();
        // Point              p3 = h->next()->next()->vertex()->point();
        // Vector             n  = CGAL::cross_product(p2-p1, p3-p1);
        CGALVector             n  = p2 - p1;
        return n / std::sqrt(n*n);
    }
};

struct movePoint
{
    CGALVector n;
    float x, y, z;

    movePoint (CGALVector normal, float dx) {
        n = normal;
        x = dx;
        y = dx;
        z = dx;
    }

    inline Point operator() (const Point& p) const {
        // TODO make it elegant
        return Point(p.x() + x*n.x(), p.y() + y*n.y(), p.z() + z*n.z());
    }
};


struct Generate_Points_at_Edges {
    const float dx_;
    FixedDistanceParticles &fdp_;

    Generate_Points_at_Edges(const float dx, FixedDistanceParticles &fdp)
        : dx_(dx), fdp_(fdp) {}

    void operator()(Facet &f) {
        // TODO use circulator
        Point A = f.halfedge()->vertex()->point();
        Point B = f.halfedge()->next()->vertex()->point();
        Point C = f.halfedge()->opposite()->vertex()->point();

        std::vector<std::pair<Point, Point>> edges {{A, B}, {A, C}, {B, C}};

        for (auto edge: edges) {
            Point F = edge.first;
            Point L = edge.second;

            float len = std::sqrt((L - F).squared_length());
            int one = 1;
            int n_points = (int) len / dx_;
            n_points = std::max(one, n_points);
            for (int j=0;j<1;j++) {
              const CGALVector N = Compute_Facet_Normal()(f);
              const CGALVector C = cross_product(L-F, N);
              for (int i = 0; i <= n_points; i++) {
                  float lambda = (float)i / ((float)n_points + 1.0);
                  // TODO use emplace back
                  // TODO shift outwards by dx in edge normal direction
                  fdp_.points.push_back(Point(F + lambda * (L - F))); // + 0.01*dx_*(float)j*C));
                  fdp_.fixId.push_back(0);
                  fdp_.maxDx.push_back(0.0);
                  fdp_.mType.push_back(0);   // Fixed Point atm
                  fdp_.dir.push_back(L - F); // Fixed Point atm
                  fdp_.facets.push_back(&f); // Fixed Point atm
              }
            }
        }
    }
};


struct Generate_Random_Points_at_Facets
{

    const float dx_;
    const float dx2_;
    const float oversubs_;
    std::vector<Point> & points_;
    // std::vector<size_t>& number_points_facet_;
    std::vector<Facet_handle> & initial_facets_;
    // std::vector<CGALVector>& facet_vectors_;

    size_t facetId;

    Generate_Random_Points_at_Facets (
            const float dx,
            const float oversubs,
            std::vector<Point> & points,
            // std::vector<size_t>& number_points_facet,
            std::vector<Facet_handle> & initial_facets
            // std::vector<CGALVector>& facet_vectors
            )
        :
        dx_(dx),
        dx2_(dx*dx),
        oversubs_(oversubs),
        points_(points),
        // number_points_facet_(number_points_facet),
        initial_facets_(initial_facets),
        // facet_vectors_(facet_vectors),
        facetId(0)
    {
        // Assume that a particle covers
        // a square area. Needed for calculation
        // of n_points
        // float dx2 = dx*dx;
        // A particle covers pi*(dx/2)^2 + (dx^2 - 2 *pi*(dx/2)^2)/2
    }

    void operator()(Facet& facet) {

        // TODO currently the Facet_iterator
        // is dereferenced and copied, this is
        // unnecessary
        // const Facet facet = Facet(*f);
        // initial_facets.push_back(facet);
        const float facet_area = Compute_Facet_Area()(facet);
        const CGALVector facet_normal =  Compute_Facet_Normal()(facet);

        // TODO make it a multiple of 4
        // for loop unrolling
        const size_t one = 1;
        // TODO compare facet area to square area

        float part_area = M_PI*dx2_/4.0;
        float add_area = dx2_/2 - part_area;
        float chance = facet_area/(part_area+add_area)*(one+oversubs_);
        int v1 = rand() % 100;
        if (chance < 1.0) {
          if ((chance  * (float) v1) > 50.0) chance = 1.0;
          else chance=0.0;
        };

        const size_t n_points = (size_t) chance;
        // std::vector<Point> points(n_points);

        Random_points_in_triangle_3<Point> g(facetToTriangle(facet));

        for (size_t i=0; i<n_points; i++) {
            // Generate a Point from iterator *g
            points_.push_back(Point(*g));
            initial_facets_.push_back(&facet);
            *g++;
        }

        // number_points_facet_[facetId] = n_points;
        // facet_vectors_[facetId] = facet_normal;

        // for (Point & elem : points) {
        //     elem = *g;
        //     *g++;
        // }

        // std::copy_n(g, n_points, std::back_inserter(paf.points)),

        // return Points_at_Facet {
        //     points,
        //     Compute_Facet_Directions()(facet),
        //     Compute_Facet_Directions()(facet),
        //     facet
        // };
        facetId++;
    }
};

struct EdgeNormal {
  CGALVector EN;
  CGALVector FN;
};


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

struct Generate_Points_at_Facets
{

    const float dx_;
    // const float dx2_;
    const float oversubs_;
    std::vector<Point> & points_;
    // std::vector<size_t>& number_points_facet_;
    std::vector<Facet_handle> & initial_facets_;
    // std::vector<CGALVector>& facet_vectors_;

    size_t facetId;

    Generate_Points_at_Facets (
            const float dx,
            const float oversubs,
            std::vector<Point> & points,
            // std::vector<size_t>& number_points_facet,
            std::vector<Facet_handle> & initial_facets
            // std::vector<CGALVector>& facet_vectors
            )
        :
        dx_(dx),
        // dx2_(dx*dx),
        oversubs_(oversubs),
        points_(points),
        // number_points_facet_(number_points_facet),
        initial_facets_(initial_facets),
        // facet_vectors_(facet_vectors),
        facetId(0)
    {
    }

    void operator()(Facet& f) {

      Point A = f.halfedge()->vertex()->point();
      Point B = f.halfedge()->next()->vertex()->point();
      Point C = f.halfedge()->opposite()->vertex()->point();
      Point S,E; // Start and end points

      const CGALVector facet_normal =  Compute_Facet_Normal()(f);
      const CGALVector edge_normal =  inwardPointingEdgeNormal(A, B, f).EN;
      const CGALVector CA = C-A;
      const CGALVector CB = C-B;

      float dist_base = dx_/2.0; // Cartesian distance from base edge
      float dist_points = dx_/2.0; // Cartesian dist from beginning of point line
      float lambda_s, lambda_e;

      float angle_s = normalise(CA)*normalise(edge_normal);
      float angle_e = normalise(CB)*normalise(edge_normal);
      size_t n_points = 0;
      int line_ctr = 0;

      do {
        lambda_s = dist_base/angle_s;
        S = A + lambda_s * normalise(CA);

        // TODO correct this
        lambda_e = dist_base/angle_e;
        E = S + lambda_e * normalise(CB);

        const CGALVector SE = S-E;

        do {

          Point P = S+normalise(SE)*dist_points;
          // std::cout
          //   << " A " << A
          //   << " B " << B
          //   << " C " << C
          //   << " NAB " << edge_normal
          //   << " S " << S
          //   << " E " << E
          //   << " dist_points " << dist_points
          //   << " P " << P
          //   << std::endl;

          dist_points += 1.5*dx_;

          n_points++;
          points_.push_back(P);
          initial_facets_.push_back(&f);
        } while (dist_points*dist_points <  (SE).squared_length()-(dx_*dx_));

        if (line_ctr % 2 == 0)  dist_points = dx_/2.0;
        else dist_points = 0;
        line_ctr++;
        dist_base += dx_/2.0;

      } while (dist_base*dist_base < (angle_s*CA).squared_length()-(dx_*dx_));

        // number_points_facet_[facetId] = n_points;
        // initial_facets_[facetId] = &f;
        // facet_vectors_[facetId] = facet_normal;

        facetId++;
    }
};

// FixedDistanceParticles create_extruded_points(
//     std::vector<Point> &orig_points,
//     int n_extrusions,
//     float dx,
//     Generate_Points_at_Facets &gpf
//                                               ) {

//     std::vector<Point> ret_p;
//     ret_p.reserve(n_extrusions * orig_points.size());

//     std::vector<size_t> ret_id;
//     ret_id.reserve(n_extrusions * orig_points.size());

//     std::vector<float> ret_dx;
//     ret_dx.reserve(n_extrusions * orig_points.size());

//     std::vector<int> ret_mType;
//     ret_mType.reserve(n_extrusions * orig_points.size());

//     std::vector<CGALVector> ret_dir;
//     ret_dir.reserve(n_extrusions * orig_points.size());

//     std::vector<Facet_handle> ret_facets;
//     ret_facets.reserve(n_extrusions * orig_points.size());

//     for (int i = 0; i < n_extrusions; i++) {
//         size_t facet_id = 0;
//         size_t points_at_facet;
//         size_t points_at_facet_ctr=0;

//         for (size_t j = 0; j < orig_points.size(); j++) {
//             // TODO iterate as long as facet handle is unchanged
//             points_at_facet = 0;

//             if (points_at_facet_ctr == points_at_facet) {
//               facet_id++;
//               points_at_facet_ctr = 0;
//             }

//             CGALVector extrude_dir = -gpf.facet_vectors_[facet_id];
//             // TODO implement shift
//             ret_id.push_back(j);
//             ret_p.push_back(
//                 Point(
//                     orig_points[j].x(),
//                     orig_points[j].y(),
//                     orig_points[j].z()) +
//                 i * 0.5*dx * extrude_dir);
//             ret_dir.push_back(extrude_dir);
//             int mType = (i==0)? 2: 3;
//             ret_mType.push_back(mType);
//             float mdx = (mType == 2)? 3.0*dx: (float)i*dx;
//             ret_dx.push_back(mdx);
//             ret_facets.push_back(gpf.initial_facets_[facet_id]);

//             points_at_facet_ctr++;
//         }
//     }

//     return {ret_p, ret_id, ret_dx, ret_mType, ret_dir, ret_facets};
// }

#endif
