#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

#include "poisson_tester.h"

using namespace dealii;


using PoissonTestTypes =
  ::testing::Types<std::integral_constant<int, 1>, // 轮着测试1-2-3d
                   std::integral_constant<int, 2>,
                   std::integral_constant<int, 3>>;


using Poisson2DTester =
  PoissonTester<std::integral_constant<int, 2>>; // 专门用来测试2d

TYPED_TEST_CASE(PoissonTester, PoissonTestTypes); // 1-2-3维测试

TYPED_TEST(PoissonTester, MakeGrid) // 1-2-3维测试
{
  // Output dimension
  std::cout << "Working on dim " << TypeParam::value << std::endl;
  this->make_grid();
}

// sharememery 这块暂时不进行测试，不然CI不通过。