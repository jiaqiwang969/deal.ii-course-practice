/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 *          Luca Heltai, 2021
 */

// Make sure we don't redefine things
#ifndef poisson_include_file
#  define poisson_include_file

#  include <deal.II/base/function.h>
#  include <deal.II/base/function_parser.h>
#  include <deal.II/base/parameter_acceptor.h>
#  include <deal.II/base/parsed_convergence_table.h>
#  include <deal.II/base/quadrature_lib.h>

#  include <deal.II/dofs/dof_handler.h>
#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/fe/fe_q.h>
#  include <deal.II/fe/fe_values.h>

#  include <deal.II/grid/grid_generator.h>
#  include <deal.II/grid/tria.h>

#  include <deal.II/lac/affine_constraints.h> // 存储所有约束, 用于处理挂起的节点
#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/precondition.h>
#  include <deal.II/lac/solver_cg.h>
#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <deal.II/numerics/data_out.h>
#  include <deal.II/numerics/matrix_tools.h>
#  include <deal.II/numerics/vector_tools.h>

#  include <fstream>
#  include <iostream>

// Forward declare the tester class
template <typename Integral>
class PoissonTester;

using namespace dealii;

// 添加模版类
template <int dim>
class Poisson : ParameterAcceptor
{
public:
  Poisson();
  void
  run();

  void
  initialize(const std::string &filename);

  void
  parse_string(const std::string &par); // 为了直接从test中嵌入代码,进行相关测试

protected:
  void
  make_grid();
  void
  refine_grid(); // 网格加密
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results(
    const unsigned cycle) const; //输出网格（不同cycle的加密方式的文件名）

  Triangulation<dim>         triangulation;
  std::unique_ptr<FE_Q<dim>> fe;
  DoFHandler<dim>            dof_handler;
  AffineConstraints<double> constraints; // 存储所有约束, 用于处理挂起的节点
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;

  FunctionParser<dim> forcing_term;                 // parameter 传呼
  FunctionParser<dim> dirichlet_boundary_condition; // 增加自定义边界条件
  FunctionParser<dim> neumann_boundary_condition; // （需要通过积分实现）


  unsigned int fe_degree           = 1;         // 默认初始参数， 1阶fe
  unsigned int n_refinements       = 4;         // 网格加密
  unsigned int n_refinement_cycles = 1;         // 误差评估再网格划分
  std::string  output_filename     = "poisson"; // 输出文件名

  std::set<types::boundary_id> dirichlet_ids = {0}; // 边界点的id标签
  std::set<types::boundary_id> neumann_ids;         // 边界点的id标签

  std::string forcing_term_expression = "1"; // 外力的表达式
  std::string dirichlet_boundary_conditions_expression = "0"; // 边界1的表达式
  std::string neumann_boundary_conditions_expression = "0"; // 边界2的表达式
  std::map<std::string, double> constants; // 恒定变量参数表

  std::string grid_generator_function  = "hyper_cube";  // 通用形状
  std::string grid_generator_arguments = "0: 1: false"; // 边界序列

  ParsedConvergenceTable error_table;

  template <typename Integral>
  friend class PossionTester; //开放测试权限
};

#endif


template class Poisson<1>; // 实例化模版类
template class Poisson<2>;
template class Poisson<3>;
