#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/write_off_points.h>

#include <iostream>
#include <fstream>
using namespace CGAL;
typedef Simple_cartesian<double>                           K;
typedef CGAL::Polyhedron_3<K>                              Polyhedron;
typedef K::Point_3                                         Point;
typedef K::FT                                              FT;

template <class PolyhedronTraits_3>
std::ofstream& operator<<( std::ofstream& out, const Polyhedron& P);

template <class PolyhedronTraits_3>
std::ifstream& operator>>( std::ifstream& in, Polyhedron& P);

int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream istream("in.off");
  istream >> polyhedron;

  int n_points = 1000000;

  // Create the generator, input is the Polyhedron polyhedron
  Random_points_in_triangle_mesh_3<Polyhedron>
      g(polyhedron);

  // Get 100 random points in cdt
  CGAL::cpp11::copy_n(g, n_points, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == n_points);

  // print the first point that was generated
  std::ofstream ostream("out.off");
  write_off_points( ostream, points);
  return 0;
}

