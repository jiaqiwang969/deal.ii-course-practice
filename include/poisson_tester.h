#ifndef poisson_tester_h
#define poisson_tester_h

#include <gtest/gtest.h>

#include <fstream>

#include "poisson.h"

using namespace dealii;

//  泊松问题的测试，使用积分常数
template <class Integral>
class PoissonTester : public ::testing::Test, public Poisson<Integral::value>
{
public:
  PoissonTester() = default;
};

#endif