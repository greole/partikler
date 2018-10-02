#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
// #include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/STL_reader.h>

#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
#include <boost/iterator/transform_iterator.hpp>

#include "git.h"
using namespace CGAL;
typedef Simple_cartesian<double>          K;
typedef CGAL::Polyhedron_3<K>             Polyhedron;
typedef K::Point_3                        Point;
typedef K::FT                             FT;
typedef Polyhedron::Facet_iterator        Facet_iterator;
typedef Polyhedron::HalfedgeDS            HalfedgeDS;

template <class PolyhedronTraits_3>
std::ofstream& operator<<( std::ofstream& out, const Polyhedron& P);

template <class PolyhedronTraits_3>
std::ifstream& operator>>( std::ifstream& in, Polyhedron& P);

struct Compute_area
{
    double operator()(const Polyhedron::Facet& f) const {
    return K::Compute_area_3()(
      f.halfedge()->vertex()->point(),
      f.halfedge()->next()->vertex()->point(),
      f.halfedge()->opposite()->vertex()->point() );
    }
};


int main(int argc, char* argv[])
{
    if(GIT_RETRIEVED_STATE) {
        std::cout << "INFO: " << GIT_HEAD_SHA1 << std::endl;
        if(GIT_IS_DIRTY) std::cerr << "WARN: there were uncommitted changes." << std::endl;
    }
    else {
        std::cerr << "WARN: failed to get the current git state. Is this a git repo?" << std::endl;
        return 0;
    }

    // Generated points are in that vector
    std::vector<Point> points;

    float dx = atof(argv[3]);

    // Assume that a particle covers
    // a square area. Needed for calculation
    // of n_points
    float dx2 = dx*dx;

    // TODO Flush cout buffer
    std::cout << "Reading input file: " <<  argv[1] << std::flush;

    std::ifstream istream(argv[1]);
    // std::vector<cpp11::array<double,3> > points,
    // std::vector<cpp11::array<int,3> > facets,
    // read_STL(istream, points, facets, false);
    Polyhedron_builder_from_STL<HalfedgeDS> builder(istream);
    std::cout << " [done]" << std::endl;

    std::cout << "Constructing polyhedron" << std::flush;
    // Create input polyhedron
    Polyhedron polyhedron;
    polyhedron.delegate(builder);
    std::cout << " [done]" << std::endl;

    std::cout << "Computing total surface area" << std::flush;
    Compute_area ca;
    float total_surface_area =
        std::accumulate(
                  boost::make_transform_iterator(polyhedron.facets_begin(), ca),
                  boost::make_transform_iterator(polyhedron.facets_end(), ca),
                  0.);
    std::cout << " [done] " << total_surface_area << std::endl;

    int n_points {total_surface_area/dx2};
    std::cout << " Number of particles " << n_points << std::endl;

    // Create the generator, input is the Polyhedron polyhedron
    Random_points_in_triangle_mesh_3<Polyhedron> g(polyhedron);

    // Get 100 random points in cdt
    CGAL::cpp11::copy_n(g, n_points, std::back_inserter(points));

    // Check that we have really created 100 points.
    assert( points.size() == n_points);

    // print the first point that was generated
    std::ofstream ostream(argv[2]);
    write_off_points( ostream, points);
    return 0;
}

