#include "FieldOps.hpp"

template <class Field, class El>
void clamp_field_in_range(Field &f, El lo, El hi) {
    for (auto &el : f) {
        el = std::clamp(el, lo, hi);
    }
}

template <> void clamp_field_in_range(FloatField &f, Scalar lo, Scalar hi) {
    for (auto &el : f) {
        el = std::clamp(el, lo, hi);
    }
}

template <> void clamp_field_in_range(VectorField &f, Scalar lo, Scalar hi) {}
