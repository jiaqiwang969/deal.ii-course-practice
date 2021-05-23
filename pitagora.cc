#include <deal.II/base/point.h>

#include <gtest/gtest.h>

using namespace dealii;


TEST(Pitagora, Norm)
{
  Point<2> x(3, 4);
  ASSERT_EQ(x.norm(), 5);
}


TEST(Pitagora, Distance)
{
  Point<2> x(4, 5);

  Point<2> y(1, 1);
  ASSERT_EQ(x.distance(y), 5);
}


int
main(int argc, char *argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
