#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
// #include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/STL_reader.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
// #include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range.hpp>

#include "git.h"
#include "io.h"

/* file IO */
#include <sys/stat.h>
#include <pthread.h>

using namespace CGAL;
typedef Simple_cartesian<double>          K;
typedef CGAL::Polyhedron_3<K>             Polyhedron;
typedef K::Point_3                        Point;
typedef Polyhedron::Traits::Vector_3      Vector;
typedef K::FT                             FT;
typedef Polyhedron::Facet_iterator        Facet_iterator;
typedef Polyhedron::HalfedgeDS            HalfedgeDS;
typedef Polyhedron::Facet                 Facet;
typedef Polyhedron::Vertex                Vertex;
typedef Polyhedron::Halfedge_const_handle HalfedgeConstHandle;
// typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

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

struct movePoint
{
    Vector n;
    float x;

    movePoint (Vector normal, float dx) {
        n = normal;
        x = dx;
    }

    inline Point operator() (const Point& p) const {
        // TODO make it elegant
        return Point(p.x() + x*n.x(), p.y() + x*n.y(), p.z() + x*n.z());
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

    std::cout << "Reading input file: " <<  argv[1] << std::flush;

    std::ifstream istream(argv[1]);

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

     for ( Facet_iterator facet_ptr = polyhedron.facets_begin();
             facet_ptr != polyhedron.facets_end();
             ++facet_ptr) {

        const Facet facet = Facet(*facet_ptr);
        float facet_area = ca(facet);
        Vector facet_normal =  cn(facet);
        int n_points {facet_area/dx2};

        // Create the generator, input is the Polyhedron polyhedron
        // TODO: get the facet constructor to work
        // Random_points_in_triangle_3<Point> g(facet);
        std::vector<Point> tmp_points;
        Random_points_in_triangle_3<Point> g(
                facet.halfedge()->vertex()->point(),
                facet.halfedge()->next()->vertex()->point(),
                facet.halfedge()->opposite()->vertex()->point()
                );

        CGAL::cpp11::copy_n(g, n_points, std::back_inserter(tmp_points));
        // movePoint mp(facet_normal, dx*i);
        // std::transform(tmp_points.begin(), tmp_points.end(), tmp_points.begin(), mp);

        for(auto const& point: tmp_points) {points.push_back(Point(point));}
        Tree tree(
            CGAL::faces(polyhedron).first,
            CGAL::faces(polyhedron).second,
            polyhedron);

        // TODO check what it does?
        // tree.accelerate_distance_queries()

        // closest points
        std::vector<Point> closest_points;
        for(auto const& point: points) {
            Point_and_primitive_id pp = tree.closest_point_and_primitive(point);
            Polyhedron::Face_handle f = pp.second; // closest primitive id

            points.push_back(
                    tree.closest_point(point));
        }
    }

    // print the first point that was generated
    std::cout << "Writing output" << std::endl;
    std::ofstream ostream(argv[2]);
    write_off_points(ostream, points);

    /** write .sph Format */
    for (int i =0; i<10; i++){
        writeData_SPH("daten", i, points);
    }
    return 0;
}
