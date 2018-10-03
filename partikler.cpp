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
#include <boost/range.hpp>

#include "git.h"

using namespace CGAL;
typedef Simple_cartesian<double>          K;
typedef CGAL::Polyhedron_3<K>             Polyhedron;
typedef K::Point_3                        Point;
typedef Polyhedron::Traits::Vector_3      Vector;
typedef K::FT                             FT;
typedef Polyhedron::Facet_iterator        Facet_iterator;
typedef Polyhedron::HalfedgeDS            HalfedgeDS;
typedef Polyhedron::Facet                 Facet;
typedef Polyhedron::Halfedge_const_handle HalfedgeConstHandle;

template <class PolyhedronTraits_3>
std::ofstream& operator<<( std::ofstream& out, const Polyhedron& P);

template <class PolyhedronTraits_3>
std::ifstream& operator>>( std::ifstream& in, Polyhedron& P);

struct Compute_area
{
    double operator()(const Facet& f) const {
    return K::Compute_area_3()(
      f.halfedge()->vertex()->point(),
      f.halfedge()->next()->vertex()->point(),
      f.halfedge()->opposite()->vertex()->point() );
    }
};

struct ComputeFacetNormal
{
    // Compute normal of the given facet.
    // Facet can be triangle, quadrilateral or a polygon as long as its planar.
    // Use first three vertices to compute the normal.
    inline Vector operator() (const Facet& f) const {
        HalfedgeConstHandle h  = f.halfedge();
        Point              p1 = h->vertex()->point();
        Point              p2 = h->next()->vertex()->point();
        Point              p3 = h->next()->next()->vertex()->point();
        Vector             n  = CGAL::cross_product(p2-p1, p3-p1);
        return n / std::sqrt(n*n);
    }
};


int main(int argc, char* argv[]) {
    if(GIT_RETRIEVED_STATE) {
        std::cout
            << "INFO: "
            << GIT_HEAD_SHA1
            << std::endl;
        if(GIT_IS_DIRTY) {
            std::cerr
                << "WARN: there were uncommitted changes."
                << std::endl;
        }
    }
    else {
        std::cerr
            << "WARN: failed to get the current git state."
            << "Is this a git repo?"
            << std::endl;
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

    Compute_area ca;
    ComputeFacetNormal cn;
    // std::vector<double> facet_areas;
    // std::vector<Vector> facet_normals;
    // std::cout << "Computing facet areas" << std::flush;
    // std::transform(
    //         polyhedron.facets_begin(),
    //         polyhedron.facets_end(),
    //         std::back_inserter(facet_areas),
    //         ca);
    // std::cout << " [done] " << std::endl;
    //
    // std::cout << "Computing facet normals" << std::flush;
    // std::transform(
    //         polyhedron.facets_begin(),
    //         polyhedron.facets_end(),
    //         std::back_inserter(facet_normals),
    //         cn);
    // std::cout << " [done] " << std::endl;

     for ( Facet_iterator facet_ptr = polyhedron.facets_begin();
             facet_ptr != polyhedron.facets_end();
             ++facet_ptr) {

        float facet_area = ca(*facet_ptr);
        Vector facet_normal =  cn(*facet_ptr);
        int n_points {facet_area/dx2};

        // Create the generator, input is the Polyhedron polyhedron
        Random_points_in_triangle_3<Point> g(*facet_ptr);

        CGAL::cpp11::copy_n(g, n_points, std::back_inserter(points));

        // Check that we have really created 100 points.
        assert( points.size() == n_points);
    }

    // print the first point that was generated
    std::cout << "Writing output" << std::flush;
    std::ofstream ostream(argv[2]);
    write_off_points(ostream, points);

    return 0;
}

