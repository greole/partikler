#include "Field.hpp"

PointField& operator+=(PointField& a, VectorField& b){
    for(size_t i=0; i<a.size(); i++){
        a[i] += {b[i][0], b[i][1], b[i][2]};
    }
    return a;
}

// template <> void IntField::write_to_disk(std::string path) {
//     write_scalar_field_impl(*this, path, get_name(), "int", "32");
// }

// template <> void FloatField::write_to_disk(std::string path) {
//     write_scalar_field_impl(*this, path, get_name(), "float", "32");
// }

// template <> void VectorField::write_to_disk(std::string path) {
//     write_vector_field_impl(
//         *this, path, get_name(), get_comp_names(), "float", "32");
// }

// template <> void PointField::write_to_disk(std::string path) {
//     write_vector_field_impl(
//         *this, path, get_name(), get_comp_names(), "float", "32");
// }
