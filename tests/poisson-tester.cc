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


// 仅对二维测试1：
TEST_F(Poisson2DTester, TestLinear)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x" << std::endl
      << "  set Dirichlet boundary ids                  = 0" << std::endl
      << "  set Finite element degree                   = 1" << std::endl
      << "  set Forcing term expression                 = 0" << std::endl
      << "  set Grid generator arguments                = 0: 1: false"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = " << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = poisson" << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str()); // 将string转化为参数parameter
  make_grid();
  setup_system();
  assemble_system();
  solve();

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution; // 与参考结果对比，对possion方程，有解析解可以作为对比

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10); // 结果评估
}

// 二维测试2:
TEST_F(Poisson2DTester, TestQuadratic)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x^2" << std::endl
      << "  set Dirichlet boundary ids                  = 0" << std::endl
      << "  set Finite element degree                   = 2" << std::endl
      << "  set Forcing term expression                 = -2" << std::endl
      << "  set Grid generator arguments                = 0: 1: false"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = " << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = quadratic"
      << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();
  setup_system();
  assemble_system();
  solve();

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}


// 二维测试3:
TEST_F(Poisson2DTester, TestMixedBC1)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x^2" << std::endl
      << "  set Dirichlet boundary ids                  = 1,2,3" << std::endl
      << "  set Finite element degree                   = 2" << std::endl
      << "  set Forcing term expression                 = -2" << std::endl
      << "  set Grid generator arguments                = 0: 1: true"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = 0" << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = quadratic"
      << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();
  setup_system();
  assemble_system();
  solve();

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}



// 二维测试：HangingNodes 为后面的自适应网格优化做铺垫
TEST_F(Poisson2DTester, TestLinearWithHangingNodes)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x" << std::endl
      << "  set Dirichlet boundary ids                  = 0" << std::endl
      << "  set Finite element degree                   = 1" << std::endl
      << "  set Forcing term expression                 = 0" << std::endl
      << "  set Grid generator arguments                = 0: 1: false"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = " << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 2" << std::endl
      << "  set Output filename                         = lin_with_handing"
      << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();
  for (unsigned int i = 0; i < 2; ++i)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->center().square() <= .25)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }


  setup_system();
  assemble_system();
  solve();
  output_results(0);

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}
