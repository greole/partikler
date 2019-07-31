#ifndef CGALTYPEDEFS_H
#define CGALTYPEDEFS_H

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
#include <CGAL/bounding_box.h>

#include <array>

using namespace CGAL;
typedef Simple_cartesian<double>          K;
typedef CGAL::Polyhedron_3<K>             CGALPolyhedron;
typedef CGAL::Triangle_3<K>               Triangle;
typedef CGAL::Aff_transformation_3<K>     Transformation;
typedef K::Point_3                        Point;
typedef K::Plane_3                        Plane;
typedef K::Line_3                         Line;
typedef CGALPolyhedron::Traits::Vector_3      CGALVector;
typedef std::array<float, 3>              Vector;

typedef K::FT                             FT;

typedef CGALPolyhedron::Facet_iterator        Facet_iterator;
typedef CGALPolyhedron::Facet_handle          Facet_handle;
typedef CGALPolyhedron::HalfedgeDS            HalfedgeDS;
typedef CGALPolyhedron::Facet                 Facet;
typedef CGALPolyhedron::Vertex                Vertex;
typedef CGALPolyhedron::Halfedge_const_handle HalfedgeConstHandle;
// typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_face_graph_triangle_primitive<CGALPolyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

template <class PolyhedronTraits_3>
std::ofstream& operator<<( std::ofstream& out, const CGALPolyhedron& P);

template <class PolyhedronTraits_3>
std::ifstream& operator>>( std::ifstream& in, CGALPolyhedron& P);

#endif
