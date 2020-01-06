#include "cgal/CGALHelper.hpp"
#include "cgal/CGALTYPEDEFS.hpp"
#include "core.hpp"
#include "gtest/gtest.h"
#include <vector>

TEST(FloatField, reoderVector) {

    std::vector<size_t> idx {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    FloatField v(10, 0.0);

    for (size_t i = 0; i < idx.size(); i++) {
        v[i] = float(i);
    }

    reorder_vector(v, idx);

    for (size_t i = 0; i < idx.size(); i++) {
        ASSERT_EQ(v[i], float(9 - i));
    }
}

TEST(FloatField, copyStdVector) {

    std::vector<float> a(2, 2.0);

    FloatField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);
}

TEST(FloatField, eagerCopy) {

    FloatField a(2, 2.0);

    FloatField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);
}

TEST(FloatField, sumABTest) {

    FloatField a(2, 1.0);

    FloatField b(2, 2.0);

    NeighbourFieldAB n(1, {0, 1});

    FloatField res(2, 0.0);

    sum_AB_impl(res, n, A(a) + B(b));

    ASSERT_EQ(res[0], 3.0);
    ASSERT_EQ(res[1], 3.0);
}

TEST(FloatField, ComplexArithmeticFields) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, (a + b) / b);

    ASSERT_EQ(res[0], 1.5);
}

TEST(FloatField, MultFloatFields) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, a * b);

    ASSERT_EQ(res[0], 2.0);
}

TEST(FloatField, AddFloatFields) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, a + b);

    ASSERT_EQ(res[0], 3.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
