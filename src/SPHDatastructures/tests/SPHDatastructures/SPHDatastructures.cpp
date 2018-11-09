#include "gtest/gtest.h"
#include "SPHDatastructures.h"
#include "CGALTYPEDEFS.h"
#include <vector>

TEST(SPHDatastructures, InplaceAdd) {

    std::vector<float> a_f { 1., 2., 3., 4.};
    std::vector<float> b_f { 1., 2., 3., 4.};
    SPHScalarField a {a_f};
    SPHScalarField b {b_f};

    a+=b;

    for (size_t ctr=0; ctr<4; ctr++) {
        std::cout << a.get_field()[ctr] << std::endl;
    }
}

TEST(SPHDatastructures, eAdd) {

    std::vector<float> a_f { 1., 2., 3., 4.};
    std::vector<float> b_f { 1., 2., 3., 4.};
    SPHScalarField a {a_f};
    SPHScalarField b {b_f};

    SPHScalarField c = a+b;

    for (size_t ctr=0; ctr<4; ctr++) {
        std::cout << c.get_field()[ctr] << std::endl;
    }
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
