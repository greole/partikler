#include <CGALTYPEDEFS.h>
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

        // TODO make it a multiple of 4
        // for loop unrolling
        const size_t n_points {(size_t) (facet_area/dx2_)};
        // std::vector<Point> points(n_points);

        Random_points_in_triangle_3<Point> g(facetToTriangle(facet));

        for (size_t i=0; i<n_points; i++) {
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
