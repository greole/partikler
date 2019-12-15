#include "cgal/CGALTYPEDEFS.hpp"
#include "gtest/gtest.h"
#include "cgal/CGALHelper.hpp"
#include "core.hpp"
#include <vector>

TEST (FloatField, copyStdVector ) {

    std::vector<float> a(2, 2.0);

    FloatField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);

}

TEST (FloatField, eagerCopy ) {

    FloatField a(2, 2.0);

    FloatField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);

}


TEST (FloatField, sumABTest ) {

    FloatField a(2, 1.0);

    FloatField b(2, 2.0);

    NeighbourFieldAB n(1, {0, 1});

    FloatField res(2, 0.0);

    sum_AB(res, n, A(a) + B(b));

    ASSERT_EQ(res[0], 3.0);
    ASSERT_EQ(res[1], 3.0);

}

TEST (FloatField, ComplexArithmeticFields ) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, (a+b)/b);

    ASSERT_EQ(res[0], 1.5);

}


TEST (FloatField, MultFloatFields ) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, a*b);

    ASSERT_EQ(res[0], 2.0);

}

TEST (FloatField, AddFloatFields ) {

  FloatField a(1, 1.0);

  FloatField b(1, 2.0);

  FloatField res(1, 0.0);

  solve(res, a+b);

  ASSERT_EQ(res[0], 3.0);

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
