#ifndef linear_elasticity_tester_h
#define linear_elasticity_tester_h

#include <gtest/gtest.h>

#include <fstream>

#include "linear_elasticity.h"

using namespace dealii;

// Test Fixture for LinearElasticity problem, using integralconstant
template <class Integral>
class LinearElasticityTester : public ::testing::Test,
                               public LinearElasticity<Integral::value>
{
public:
  LinearElasticityTester() = default;
};

#endif