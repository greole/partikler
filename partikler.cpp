#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range.hpp>

// #include <future>
// #include <execution>
#include <parallel/algorithm>
#include <parallel/settings.h>
#include <omp.h>

#include "git.h"
#include "io.h"
#include "particle_helper.h"
#include "SearchCubes.h"

/* file IO */
#include <sys/stat.h>
#include <pthread.h>
#include <CGAL_TYPEDEFS.h>

int main(int argc, char* argv[]) {

    print_git_state();

    __gnu_parallel::_Settings s;
    s.algorithm_strategy = __gnu_parallel::force_sequential;
    __gnu_parallel::_Settings::set(s);

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

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    const float dx = atof(argv[3]);

    // Explicitly set number of threads.
    // const int threads_wanted = 4;
    // omp_set_dynamic(false);
    // omp_set_num_threads(threads_wanted);


    // s.algorithm_strategy = __gnu_parallel::force_parallel;
    // __gnu_parallel::_Settings::set(s);

    const size_t number_of_facets = std::distance(
            polyhedron.facets_begin(),
            polyhedron.facets_end());

    std::cout << "number_of_facets" << number_of_facets << std::endl;
    //TODO std::vector:reserve
    std::vector<Point>   points;
    std::vector<size_t>  number_points_facet(number_of_facets);
    std::vector<Facet>   initial_facets(number_of_facets);
    std::vector<Vector>  facet_vectors(number_of_facets);

    Generate_Points_at_Facets gpf(
            dx, points,
            number_points_facet,
            initial_facets,
            facet_vectors
    );

    std::for_each(
        polyhedron.facets_begin(),
        polyhedron.facets_end(),
        gpf
    );

    const size_t n_points = points.size();

    std::cout << "number_of_points:" <<  n_points << std::endl;

    // TODO reimplement layer extrusion

    // SPH steps
    // compute distance to nearest neighbour points
    // generates 2d array distances [npoints][default_neigbours]

    // initialise search cubes
    // an array of particle ids and cgal circulator like iterator
    // search_cube[id][particle_id]
    // supports search_cube.cube(id).begin();
    // search_cubes
    //
    SearchCubeTree search_cube_tree (points, 3*dx);

    const size_t default_neigbours = 40;

    // intializes to 0

    // create_neighbours(n_points, dx, default_neigbours, search_cubes, neighbours);

    // // print the first point that was generated
    // std::cout << "Writing output" << std::endl;
    //
    // Tree tree(
    //     CGAL::faces(polyhedron).first,
    //     CGAL::faces(polyhedron).second,
    //     polyhedron);
    //
    // // TODO check what it does?
    // tree.accelerate_distance_queries();
    //
    // for (int i =0; i<10; i++){
    //
    //     // closest points
    //     std::vector<Point> transformed_points {};
    //
    //     std::cout << "Timestep " << i << std::endl;
    //     // for(auto const& point: points) {
    //     for(int id=0; id < points.size(); id++){
    //         Point point = points[id];
    //         Facet initial_facet = initial_facets[i];
    //         const Vector facet_vector = facet_vectors[i];
    //         // std::cout << "id" << id << std::endl;
    //
    //
    //         // Test if point is still in initial triangle
    //         // if not find new Triangle
    //         // TODO replace closest_point_and_primitive by
    //         // search over neighbour triangles
    //         if (! isInsideTriangle(initial_facets[id])(point)) {
    //             // std::cout
    //             //     << "Point moved out of initial facet"
    //             //     << std::endl;
    //             Point_and_primitive_id pp =
    //                 tree.closest_point_and_primitive(point);
    //             initial_facets[i] = Facet(*pp.second); // closest primitive id
    //             initial_facet = Facet(*pp.second);
    //         }
    //
    //         // TODO dont recompute normal vector
    //         // use normal vector by face_id
    //         movePoint mp(facet_vector, dx);
    //         transformed_points.push_back(mp(point));
    //     }
    //
    //     std::ofstream ostream(argv[2] + intToStr(i) + ".off");
    //     write_off_points(ostream, transformed_points);
    //     points = transformed_points;
    //     #<{(|* write .sph Format |)}>#
    //     writeData_SPH("daten", i, points);
    // }
    return 0;
}
