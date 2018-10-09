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
#include <CGAL/Triangle_3.h>

using namespace CGAL;
typedef Simple_cartesian<double>          K;
typedef CGAL::Polyhedron_3<K>             Polyhedron;
typedef CGAL::Triangle_3<K>               Triangle;
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

struct ComputeFacetVector
{
    // Compute normal of the given facet.
    // Facet can be triangle, quadrilateral or a polygon as long as its planar.
    // Use first three vertices to compute the normal.
    inline Vector operator() (const Facet& f) const {
        HalfedgeConstHandle h  = f.halfedge();
        Point              p1 = h->vertex()->point();
        Point              p2 = h->next()->vertex()->point();
        // Point              p3 = h->next()->next()->vertex()->point();
        // Vector             n  = CGAL::cross_product(p2-p1, p3-p1);
        Vector             n  = p2 - p1;
        return n / std::sqrt(n*n);
    }
};

struct movePoint
{
    Vector n;
    float x, y, z;

    movePoint (Vector normal, float dx) {
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

Triangle facetToTriangle(Polyhedron::Face_handle f) {
    const Facet facet {*f};
    return Triangle(
           facet.halfedge()->vertex()->point(),
           facet.halfedge()->next()->vertex()->point(),
           facet.halfedge()->opposite()->vertex()->point()
       );
}

struct isInsideTriangle
{
    Polyhedron::Face_handle f;

    isInsideTriangle (Polyhedron::Face_handle f) {
        f = f;
    }

    inline bool operator() (const Point& p) const {
        // NOTE implementation of triangle.has_on(p)
        return facetToTriangle(f).has_on(p);
    }

};


