#include <gtest/gtest.h>

#include <fstream>

#include "poisson.h"

using namespace dealii;

// Test Fixture for step-3
class PoissonTester : public ::testing::Test, public Poisson
{
public:
  PoissonTester() = default;
};


TEST_F(PoissonTester, MakeGrid)
{
  make_grid();
}