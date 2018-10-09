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
#include "particle_helper.h"

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
    std::vector<Polyhedron::Face_handle> initial_facet;

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

        initial_facet.push_back(facet_ptr);
        const Facet facet = Facet(*facet_ptr);
        float facet_area = ca(facet);
        Vector facet_normal =  cn(facet);
        int n_points {facet_area/dx2};

        // Create the generator, input is the Polyhedron polyhedron
        // TODO: get the facet constructor to work
        // Random_points_in_triangle_3<Point> g(facet);
        std::vector<Point> tmp_points;
        Random_points_in_triangle_3<Point> g(facetToTriangle(facet_ptr));

        std::copy_n(g, n_points, std::back_inserter(tmp_points));

        for(auto const& point: tmp_points) {points.push_back(Point(point));}
    }

    // print the first point that was generated
    std::cout << "Writing output" << std::endl;

    Tree tree(
        CGAL::faces(polyhedron).first,
        CGAL::faces(polyhedron).second,
        polyhedron);

    // TODO check what it does?
    tree.accelerate_distance_queries();

    for (int i =0; i<10; i++){

        // closest points

        std::vector<Point> transformed_points {};

        std::cout << "Timestep " << i << std::endl;
        for(auto const& point: points) {
            Point_and_primitive_id pp = tree.closest_point_and_primitive(point);
            Polyhedron::Face_handle f = pp.second; // closest primitive id
            ComputeFacetVector cm;
            Vector normal = cm(*f);
            movePoint mp(normal, dx);

            transformed_points.push_back(mp(point));
            std::ofstream ostream(argv[2] + intToStr(i) + ".off");
            write_off_points(ostream, transformed_points);
        }

        points = transformed_points;
        /** write .sph Format */
        writeData_SPH("daten", i, points);
    }
    return 0;
}
