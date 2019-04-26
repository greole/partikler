#include <CGALTYPEDEFS.hpp>
#include <math.h>

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
    inline CGALVector operator() (const Facet& f) const {
        HalfedgeConstHandle h  = f.halfedge();
        Point              p1 = h->vertex()->point();
        Point              p2 = h->next()->vertex()->point();
        Point              p3 = h->next()->next()->vertex()->point();
        CGALVector             n  = CGAL::cross_product(p2-p1, p3-p1);
        return n / std::sqrt(n*n);
    }
};

struct Compute_Facet_Directions
{
    // Compute normal of the given facet.
    // Facet can be triangle, quadrilateral or a polygon as long as its planar.
    // Use first three vertices to compute the normal.
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

        for (auto edge : edges) {
            Point F = edge.first;
            Point L = edge.second;

            float len = std::sqrt((L - F).squared_length());
            int one = 1;
            int n_points = (int)len / dx_;
            n_points = std::max(one, n_points);
            for (int j=0;j<3;j++) {
              const CGALVector N = Compute_Facet_Normal()(f);
              const CGALVector C = cross_product(L-F, N);
              for (int i = 0; i <= n_points; i++) {
                  float lambda = (float)i / ((float)n_points + 1.0);
                  // TODO use emplace back
                  fdp_.points.push_back(Point(F + lambda * (L - F) + 0.01*dx_*(float)j*C));
                  fdp_.fixId.push_back(0);
                  fdp_.maxDx.push_back(dx_);
                  fdp_.mType.push_back(0);   // Fixed Point atm
                  fdp_.dir.push_back(L - F); // Fixed Point atm
              }
            }
        }
    }
};

struct Generate_Points_at_Facets
{

    const float dx_;
    const float dx2_;
    std::vector<Point> & points_;
    std::vector<size_t>& number_points_facet_;
    std::vector<Facet_handle> & initial_facets_;
    std::vector<CGALVector>& facet_vectors_;

    size_t facetId;

    Generate_Points_at_Facets (
            const float dx,
            std::vector<Point> & points,
            std::vector<size_t>& number_points_facet,
            std::vector<Facet_handle> & initial_facets,
            std::vector<CGALVector>& facet_vectors
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
        const CGALVector facet_normal =  Compute_Facet_Normal()(facet);

        // TODO make it a multiple of 4
        // for loop unrolling
        const size_t one = 1;
        const size_t n_points = std::max(one, {(size_t) (facet_area/dx2_)});
        // std::vector<Point> points(n_points);

        Random_points_in_triangle_3<Point> g(facetToTriangle(facet));

        for (size_t i=0; i<n_points; i++) {
            // Generate a Point from iterator *g
            points_.push_back(Point(*g));
            *g++;
        }

        number_points_facet_[facetId] = n_points;
        initial_facets_[facetId] = &facet;
        facet_vectors_[facetId] = facet_normal;

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

FixedDistanceParticles create_extruded_points(
    std::vector<Point> &orig_points,
    int n_extrusions,
    float dx,
    Generate_Points_at_Facets &gpf
                                              ) {

    std::vector<Point> ret_p;
    ret_p.reserve(n_extrusions * orig_points.size());

    std::vector<size_t> ret_id;
    ret_id.reserve(n_extrusions * orig_points.size());

    std::vector<float> ret_dx;
    ret_dx.reserve(n_extrusions * orig_points.size());

    std::vector<int> ret_mType;
    ret_mType.reserve(n_extrusions * orig_points.size());

    std::vector<CGALVector> ret_dir;
    ret_dir.reserve(n_extrusions * orig_points.size());

    std::vector<Facet_handle> ret_facets;
    ret_facets.reserve(n_extrusions * orig_points.size());

    for (int i = 0; i < n_extrusions; i++) {
        size_t facet_id = 0;
        size_t points_at_facet;
        size_t points_at_facet_ctr=0;

        for (size_t j = 0; j < orig_points.size(); j++) {
            points_at_facet = gpf.number_points_facet_[facet_id];

            if (points_at_facet_ctr == points_at_facet) {
              facet_id++;
              points_at_facet_ctr = 0;
            }

            CGALVector extrude_dir = -gpf.facet_vectors_[facet_id];
            // TODO implement shift
            ret_id.push_back(j);
            ret_p.push_back(
                Point(
                    orig_points[j].x(),
                    orig_points[j].y(),
                    orig_points[j].z()) +
                i * 0.5*dx * extrude_dir);
            ret_dir.push_back(extrude_dir);
            int mType = (i==0)? 2: 3;
            ret_mType.push_back(mType);
            float mdx = (mType == 2)? 3.0*dx: (float)i*dx;
            ret_dx.push_back(mdx);
            ret_facets.push_back(gpf.initial_facets_[facet_id]);

            points_at_facet_ctr++;
        }
    }

    return {ret_p, ret_id, ret_dx, ret_mType, ret_dir, ret_facets};
}


// returns ids of all particles with at least one neighbour
std::vector<size_t> find_isolated_particles(const SortedNeighbours &neighbours){


}
