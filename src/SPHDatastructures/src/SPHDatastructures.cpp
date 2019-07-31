#include "SPHDatastructures.hpp"

STLSurfaceDist compute_STLSurface_dist(
    Point opos, Point npos,
    Facet_handle start, Facet_handle end) {

    if (start == end) {
        // if particles on same facet return Cartesian distance
        CGALVector lenVo = npos - opos;
        CGALVector lenVn = opos - npos;
        return {(float) length(lenVo), lenVo, lenVn};
    }

    std::pair<bool, Path> ret_path = searchPath({start, end});

    if (!ret_path.first) {
        // if path search failed use Cartesian distance 
        CGALVector lenVo = npos - opos;
        CGALVector lenVn = opos - npos;
        return {(float) length(lenVo), lenVo, lenVn};
    }

    Path path = ret_path.second;

    Point nposp, oposp;
    nposp = projectedPoint(path, npos);
    reverse(path.begin(), path.end());
    oposp = projectedPoint(path, opos);

    CGALVector lenVo = nposp - opos;
    CGALVector lenVn = oposp - npos;

    return {(float) length(lenVo), lenVo, lenVn};
}

SearchCubeDomain
initSearchCubeDomain(const std::vector<Point> particles, float dx)
{

    size_t n_particles = particles.size();
    auto bound_box = bounding_box(particles.begin(), particles.end());

    // bounds are scaled a bit to avoid particles on exact boundary 
    const float domain_extrusion = 0.01*dx; // in dx
    float min_x = bound_box.min().x() - domain_extrusion;
    float min_y = bound_box.min().y() - domain_extrusion;
    float min_z = bound_box.min().z() - domain_extrusion;
    float max_x = bound_box.max().x() + domain_extrusion;
    float max_y = bound_box.max().y() + domain_extrusion;
    float max_z = bound_box.max().z() + domain_extrusion;
    // std::cout << "initSearchCubeDomain" << " min_z " << min_z << std::endl;
    // std::cout << "initSearchCubeDomain" << " max_z " << max_z << std::endl;
    // for (auto p: particles) std::cout << p << std::endl;

    // ceil with fixed dx increases max(x,y,z) a bit
    unsigned int n_cubes_x = (unsigned int)ceil(((max_x - min_x) / dx));
    unsigned int n_cubes_y = (unsigned int)ceil(((max_y - min_y) / dx));
    unsigned int n_cubes_z = (unsigned int)ceil(((max_z - min_z) / dx));

    // Check if domain is not degenerated
    // search cubes currently assume 26 neighbour cubes
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_x >= 3)
        << " n_cubes_x less than 3";
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_y >= 3)
        << " n_cubes_y less than 3";
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_z >= 3)
        << " n_cubes_z less than 3";

    std::cout
      << " [DEBUG] SPHDatastructures.hpp::initSearchCubeDomain"
      << " min_x " << min_x
      << " min_y " << min_y
      << " min_z " << min_z
      << " max_x " << max_x
      << " max_y " << max_y
      << " max_z " << max_z
      << " n_cubes_x " << n_cubes_x
      << " n_cubes_y " << n_cubes_y
      << " n_cubes_z " << n_cubes_z
      << " dx " << dx
      << std::endl;

    const float idx = 1.0 / dx;
    return SearchCubeDomain {
        Point3D {min_x, min_y, min_z},
        Point3D {max_x, max_y, max_z},
        dx,
        idx,
        SubDivision {n_cubes_x, n_cubes_y, n_cubes_z},
        (size_t)n_cubes_x * (size_t)n_cubes_y * (size_t)n_cubes_z,
    };
}

template <class T>
void reorder_vector(const std::vector<size_t>& idxs, std::vector<T>& vec) {
    std::vector<T> tmp(vec.size());

    for(size_t i=0; i<idxs.size(); i++) {
        tmp[idxs[i]] = vec[i];
    }

    vec=tmp;
}


// template void reorder_vector<kFacet_handle>;
template void reorder_vector(const std::vector<size_t>&, std::vector<Facet_handle>& );
template void reorder_vector(const std::vector<size_t>&, std::vector<size_t>& );
template void reorder_vector(const std::vector<size_t>&, std::vector<NeighbourPair>& );
template void reorder_vector(const std::vector<size_t>&, std::vector<VectorPair>& );
template void reorder_vector(const std::vector<size_t>&, std::vector<STLSurfaceDist>& );
template void reorder_vector(const std::vector<size_t>&, std::vector<SearchCube>& );


SPHVectorField
particle_distance_vec(const SPHPointField &pos, const SPHField<NeighbourPair> &np) {

    const size_t ret_size {np.size()};

    std::vector<Vector> ret(ret_size);

    OP_NEIGHBOUR_LOOP(const Point &opos = pos[oid];
                      const Point &npos = pos[nid];
                      const CGALVector lenV = npos - opos;

                      ret[ctr][0] = lenV[0];
                      ret[ctr][1] = lenV[1];
                      ret[ctr][2] = lenV[2];)

        return SPHVectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true);
}

std::ostream &operator<<(std::ostream &os, SPHFloatField const &m) {

    auto fs = m.get_field();
    os
        << "\n"
        << m.get_name() << " " << fs.size()
        << std::endl;

    size_t id = 0;
    for (auto f : fs) {
        os << m.get_name() << "particle " << id << " " << f << "\n" << std::flush;
        id++;
    }
    return os;

}

std::ostream &operator<<(std::ostream &os, SPHVectorField const &m) {

    auto fs = m.get_field();
    os 
        << "\n"
        << m.get_name() << " " << fs.size()
        << std::endl;
    size_t id = 0;
    for (auto f : fs) {os
            << m.get_name()
            << " Particle " << id << " ("
            << f[0] << " "
            << f[1] << " "
            << f[2] << ")\n" << std::flush;
        id++;
    }
    return os;
}


std::ostream &operator<<(std::ostream &os, SPHPointField const &m) {

    auto fs = m.get_field();
    os
        << "\n"
        << m.get_name() << " " << fs.size()
        << std::endl;
    size_t id = 0;


    for (auto f : fs) {os
            << m.get_name()
            << " Particle " << id << " ("
            << f[0] << " "
            << f[1] << " "
            << f[2] << ")\n" << std::flush;
        id++;
    }
    return os;
}

void SPHIntField::write_to_disk(std::string path) {
    write_int_field(
        path,
        get_name(),
        f_,
        get_type(),
        "32"
        );
};

void SPHSizeTField::write_to_disk(std::string path) {
    write_int_field(
        path,
        get_name(),
        f_,
        "int",// get_type(),
        "32"//"64"
        );

};

void SPHFloatField::write_to_disk(std::string path) {
    write_field(
        path,
        get_name(),
        f_,
        get_type(),
        "32"
        );
};

void SPHVectorField::write_to_disk(std::string path) {
    write_vector_field(
        path,
        get_name(),
        f_,
        get_type(),
        "32",
        get_comp_names()
        );
};

void SPHPointField::write_to_disk(std::string path) {
    write_point_field(
        path,
        get_name(),
        f_,
        get_type(),
        "32",
        get_comp_names()
        );
        };
