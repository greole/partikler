/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

#include "Datastructures.hpp"

// template <class T>
// void reorder_vector(const std::vector<size_t>& idxs, std::vector<T>& vec) {
//     std::vector<T> tmp(vec.size());

//     for(size_t i=0; i<idxs.size(); i++) {
//         tmp[idxs[i]] = vec[i];
//     }

//     vec=tmp;
// }


// // template void reorder_vector<kFacet_handle>;
// template void reorder_vector(const std::vector<size_t>&, std::vector<Facet_handle>& );
// template void reorder_vector(const std::vector<size_t>&, std::vector<size_t>& );
// template void reorder_vector(const std::vector<size_t>&, std::vector<searchcubes::NeighbourPair>& );
// template void reorder_vector(const std::vector<size_t>&, std::vector<VectorPair>& );
// template void reorder_vector(const std::vector<size_t>&, std::vector<STLSurfaceDist>& );
// template void reorder_vector(const std::vector<size_t>&, std::vector<searchcubes::SearchCube>& );


// // VectorField
// // particle_distance_vec(const PointField &pos, const Field<searchcubes::NeighbourPair> &np) {

// //     const size_t ret_size {np.size()};

// //     std::vector<Vector> ret(ret_size);

// //     OP_NEIGHBOUR_LOOP(const Point &opos = pos[oid];
// //                       const Point &npos = pos[nid];
// //                       const CGALVector lenV = npos - opos;

// //                       ret[ctr][0] = lenV[0];
// //                       ret[ctr][1] = lenV[1];
// //                       ret[ctr][2] = lenV[2];)

// //         return VectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true);
// // }

// // std::ostream &operator<<(std::ostream &os, FloatField const &m) {

// //     auto fs = m.get_field();
// //     os
// //         << "\n"
// //         << m.get_name() << " " << fs.size()
// //         << std::endl;

// //     size_t id = 0;
// //     for (auto f : fs) {
// //         os << m.get_name() << "particle " << id << " " << f << "\n" << std::flush;
// //         id++;
// //     }
// //     return os;

// // }

// // std::ostream &operator<<(std::ostream &os, VectorField const &m) {

// //     auto fs = m.get_field();
// //     os 
// //         << "\n"
// //         << m.get_name() << " " << fs.size()
// //         << std::endl;
// //     size_t id = 0;
// //     for (auto f : fs) {os
// //             << m.get_name()
// //             << " Particle " << id << " ("
// //             << f[0] << " "
// //             << f[1] << " "
// //             << f[2] << ")\n" << std::flush;
// //         id++;
// //     }
// //     return os;
// // }




// // std::ostream &operator<<(std::ostream &os, PointField const &m) {

// //     auto fs = m.get_field();
// //     os
// //         << "\n"
// //         << m.get_name() << " " << fs.size()
// //         << std::endl;
// //     size_t id = 0;


// //     for (auto f : fs) {os
// //             << m.get_name()
// //             << " Particle " << id << " ("
// //             << f[0] << " "
// //             << f[1] << " "
// //             << f[2] << ")\n" << std::flush;
// //         id++;
// //     }
// //     return os;
// // }

// // void IntField::write_to_disk(std::string path) {
// //     write_int_field(
// //         path,
// //         get_name(),
// //         f_,
// //         get_type(),
// //         "32"
// //         );
// // };

// // void SizeTField::write_to_disk(std::string path) {
// //     write_int_field(
// //         path,
// //         get_name(),
// //         f_,
// //         "int",// get_type(),
// //         "32"//"64"
// //         );

// // };

// // void FloatField::write_to_disk(std::string path) {
// //     write_field(
// //         path,
// //         get_name(),
// //         f_,
// //         get_type(),
// //         "32"
// //         );
// // };

// // void VectorField::write_to_disk(std::string path) {
// //     write_vector_field(
// //         path,
// //         get_name(),
// //         f_,
// //         get_type(),
// //         "32",
// //         get_comp_names()
// //         );
// // };

// // void PointField::write_to_disk(std::string path) {
// //     write_point_field(
// //         path,
// //         get_name(),
// //         f_,
// //         get_type(),
// //         "32",
// //         get_comp_names()
// //         );
// //         };
