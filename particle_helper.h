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
#include <CGAL/bounding_box.h>

#include <math.h>

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

// Data Structures

// struct Points_at_Facet
// {
//     const std::vector<Point>& points;
//     const Vector face_direction_1;
//     const Vector face_direction_2;
//     const Facet& facet;
// };

Triangle facetToTriangle(const Facet& facet) {
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

struct Compute_Facet_Directions
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


struct isInsideTriangle
{
    Facet facet;

    isInsideTriangle (const Facet& f) {
        facet = f;
    }

    inline bool operator() (const Point& p) const {
        // std::cout << "Triangle test" << std::endl;
        // NOTE implementation of triangle.has_on(p)
        // TODO segfaults
        return facetToTriangle(facet).has_on(p);
        // return false;
    }

};


struct Generate_Points_at_Facets
{

    const float dx_;
    const float dx2_;
    std::vector<Point> & points_;
    std::vector<size_t>& number_points_facet_;
    std::vector<Facet> & initial_facets_;
    std::vector<Vector>& facet_vectors_;

    size_t facetId;

    Generate_Points_at_Facets (
            const float dx,
            std::vector<Point> & points,
            std::vector<size_t>& number_points_facet,
            std::vector<Facet> & initial_facets,
            std::vector<Vector>& facet_vectors
            )
        :
        dx_(dx),
        dx2_(dx*dx),
        points_(points),
        number_points_facet_(number_points_facet),
        initial_facets_(initial_facets),
        facet_vectors_(facet_vectors),
        facetId(0)
    {
        // Assume that a particle covers
        // a square area. Needed for calculation
        // of n_points
        // float dx2 = dx*dx;
    }

    void operator()(Facet& facet) {

        // TODO currently the Facet_iterator
        // is dereferenced and copied, this is
        // unnecessary
        // const Facet facet = Facet(*f);
        // initial_facets.push_back(facet);
        const float facet_area = Compute_Facet_Area()(facet);
        const Vector facet_normal =  Compute_Facet_Normal()(facet);
        const int n_points {facet_area/dx2_};
        // std::vector<Point> points(n_points);

        Random_points_in_triangle_3<Point> g(facetToTriangle(facet));

        for (int i=0; i<n_points; i++) {
            // Generate a Point from iterator *g
            points_.push_back(Point(*g));
            *g++;
            number_points_facet_[facetId] = n_points;
            initial_facets_[facetId] = facet;
            facet_vectors_[facetId] = facet_normal;
        }

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

class SearchCubeTree {

    private:

        // domain geometry
        const K::Iso_cuboid_3           bound_box_;
        const float                     dx_;
        const std::vector<size_t>       n_cubes_;

        std::vector<std::vector<size_t>> particles_in_cube;

        std::vector<size_t> original_cube_ids;

        std::vector<std::vector<bool>> non_empty_cube_neighbours;

    public:

        SearchCubeTree(const std::vector<Point>& points, const float dx):
            bound_box_(bounding_box(points.begin(), points.end())),
            dx_(dx),
            n_cubes_({
                ceil((bound_box_.max().x() - bound_box_.min().x())/dx),
                ceil((bound_box_.max().y() - bound_box_.min().y())/dx),
                ceil((bound_box_.max().z() - bound_box_.min().z())/dx)
                })
    {
        size_t tot_n_cubes = n_cubes_[0] * n_cubes_[1] * n_cubes_[2];
        std::cout
            << "Initialising the SearchCubeTree for "
            << n_cubes_[0] << "x" << n_cubes_[1] << "x" << n_cubes_[2]
            << "=" << tot_n_cubes
            <<  " search cubes" << std::endl;

        std::vector<size_t> def;
        std::vector<std::vector<size_t>> tmp_particles_in_cube (tot_n_cubes, def);
        // particles_in_cube = std::vector<std::vector<int>> (tot_n_cubes, def);

        // std::vector<size_t> particles_in_cubes_unsorted(points.size());
        for (int i=0; i<points.size(); i++) {
            const size_t cube_id = position_to_cube_id(points[i]);
            tmp_particles_in_cube[cube_id].push_back(i);
        }

        std::cout << "Removing empty search cubes" << std::endl;
        size_t cube_ctr = 0;
        for (auto pic: tmp_particles_in_cube) {
            if (pic.size() > 0) {
                particles_in_cube.push_back(pic);
                original_cube_ids.push_back(cube_ctr);
                cube_ctr++;
            }
        }
        const size_t n_non_empty_search_cubes = particles_in_cube.size();
        std::cout << n_non_empty_search_cubes << std::endl;

        std::cout << "Computing search cubes neighbours" << std::endl;

        non_empty_cube_neighbours = std::vector<std::vector<bool>>(n_non_empty_search_cubes);

        cube_ctr = 0;
        for (auto id: original_cube_ids) {
            non_empty_cube_neighbours[cube_ctr] = computeNeighbourId(original_cube_ids, id);
        }
    };

    const std::vector<bool> computeNeighbourId (
            const std::vector<size_t> & cube_ids,
            const size_t id) const {

            //  n_cubes_[1]=3
            //  Bottom        Middle        Top
            //  ^ y           ^ y           ^ y
            //  |06 07 08     |15 16 17     |24 25 26 // back
            //  |03 04 05     |12 13 14     |21 22 23 // middle
            //  |00 01 02     |09 10 11     |18 19 20 // front
            //  |-------> x   |-------> x   |-------> x
            //  left  right                           n_cubes_[0]=3
            //
            // search cube ids start from min(x) to max(z)
            // in a rhs coordinate system
            // - the right neighbour (x-axis) is +1
            // - the back neighbour  (y-axis) is +(n_cubes_[0])
            // - the upper neighbour (z-axis) is +(n_cubes_[0]*n_cubes[1])


            std::vector<int> offsets (26);

            int ctr = 0;
            for (int k: {-1, 0, +1}) {
                for (int j: {-1, 0, +1}) {
                    for (int i: {-1, 0, +1}) {
                        if ( ctr != 13) { // do nothing for the centre cube
                            offsets[ctr] = k * n_cubes_[0] * n_cubes_[1] + j * n_cubes_[0] + i;
                            ctr++;
                        }
                    }
                }
            }

            std::vector<bool> neighbours (26);

            ctr = 0;
            for (auto offset: offsets) {
                neighbours[ctr] =
                    std::find(cube_ids.begin(), cube_ids.end(), offset) != cube_ids.end();
            }

            return neighbours;
    }

    const std::vector<Point> get_neighbour_candidates(const Point& point) const {

            const size_t search_cube_id = position_to_cube_id(point);

            // FIXME
            // std::vector<Point> candidates(40);
            //
            // // Step 1. append all points in search_cube_id
            //
            // candidates.append(particles_in_cube[search_cube_id]);
            //
            // // Step 2. append all points neighbouring search cubes
            // neighbour_mask = non_empty_cube_neighbours[search_cube_id];

    }

    const size_t position_to_cube_id(const Point& p) const {
        const size_t i {(p.x() - bound_box_.min().x())/dx_};
        const size_t j {(p.y() - bound_box_.min().y())/dx_};
        const size_t k {(p.z() - bound_box_.min().z())/dx_};
        return i + n_cubes_[0]*j + n_cubes_[0]*n_cubes_[1]*k;
    }
};


class ParticleNeigbourMatrix  {

    private:

        const SearchCubeTree& search_cube_tree_;

        const std::vector<Point>& points_;

        const float dx;
        // size_t neighbours[n_points][default_neigbours];

        std::vector<std::vector<size_t>> neighbours;

    public:

        ParticleNeigbourMatrix (
            const SearchCubeTree& search_cube_tree,
            const std::vector<Point>& points,
            const float dx
        ):
            search_cube_tree_(search_cube_tree),
            points_(points),
            dx(dx)
        {

        }

        const std::vector<size_t> findNeigbours(
                const std::vector<Point>& points, const size_t id) {

            // Step 1. get corresponding search cube

            // Step 2. get all particle candidates
            const Point& point = points[id];
            search_cube_tree_.get_neighbour_candidates(point);


            // Step 3. compute all distances
            // Step 4. filter distance below threshold
            // Step 5. store distances and ids

        }

};
