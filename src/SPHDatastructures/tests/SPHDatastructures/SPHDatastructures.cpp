#include "gtest/gtest.h"
#include "SPHDatastructures.hpp"
#include "CGALTYPEDEFS.hpp"
#include <vector>

TEST(SPHDatastructures, InplaceAdd) {

    const size_t n = 1000000;
    std::vector<float> a_f (n, 1.0);
    std::vector<float> b_f (n, 2.0);
    SPHFloatField a {a_f};
    SPHFloatField b {b_f};

    a+=b;

    for (size_t ctr=0; ctr<a_f.size(); ctr++) {
        ASSERT_EQ(a.get_field()[ctr], 3.);
    }
}

TEST(SPHDatastructures, InplaceMultFloatVec) {

    const size_t n = 1'000'000;
    std::vector<float> a_f (n, 2.0);
    std::vector<float> b_f (n, 3.0);
    SPHFloatField a {a_f};
    SPHFloatField b {b_f};

    a*=b;

    for (size_t ctr=0; ctr<a_f.size(); ctr++) {
        ASSERT_EQ(a.get_field()[ctr], 6.) << ctr;
    }
}

TEST(SPHDatastructures, InplaceMultFloatVecSIMD) {

    const size_t n = 1'000'000;
    std::vector<float> a_f (n, 2.0);
    std::vector<float> b_f (n, 3.0);
    SPHFloatField a {a_f};
    SPHFloatField b {b_f};

    // a.simd_inplace_mult(b);

    for (size_t ctr=0; ctr<a_f.size(); ctr++) {
        ASSERT_EQ(a.get_field()[ctr], 6.) << ctr;
    }
}


TEST(SPHDatastructures, InplaceMultFloat) {

    const size_t n = 1000000;
    std::vector<float> a_f (n, 2.0);
    SPHFloatField a {a_f};

    a*=3.;

    for (size_t ctr=0; ctr<a_f.size(); ctr++) {
        ASSERT_EQ(a.get_field()[ctr], 6.);
    }
}


// TEST(SPHDatastructures, eAdd) {
//
//     std::vector<float> a_f { 1., 2., 3., 4.};
//     std::vector<float> b_f { 1., 2., 3., 4.};
//     SPHScalarField a {a_f};
//     SPHScalarField b {b_f};
//
//     SPHScalarField c = a+b;
//
//     for (size_t ctr=0; ctr<4; ctr++) {
//         std::cout << c.get_field()[ctr] << std::endl;
//     }
// }
//

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
