#include "cgal/CGALTYPEDEFS.hpp"
#include "gtest/gtest.h"
#include "cgal/CGALHelper.hpp"
#include "core.hpp"
#include <vector>

TEST (foo, bar) {
  std::cout << normalise(CGALVector {1,2,3}) << std::endl;

  ASSERT_TRUE(true);

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
