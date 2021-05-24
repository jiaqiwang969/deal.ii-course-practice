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
#include "poisson.h"

using namespace dealii;

template <int dim> // dim模版
Poisson<dim>::Poisson()
  : dof_handler(triangulation)
{
  // 方式1：将初始化的的自定义参数赋到变量上面
  add_parameter("Finite element degree", fe_degree);
  add_parameter("Number of global refinements", n_refinements);
  add_parameter("Output filename", output_filename);
  add_parameter("Forcing term expression", forcing_term_expression);
  add_parameter("Dirichlet boundary condition expression",
                dirichlet_boundary_conditions_expression);
  add_parameter("Neumann boundary condition expression",
                neumann_boundary_conditions_expression);

  add_parameter("Dirichlet boundary ids", dirichlet_ids);
  add_parameter("Neumann boundary ids", neumann_ids);

  add_parameter("Problem constants", constants);
  add_parameter("Grid generator function", grid_generator_function);
  add_parameter("Grid generator arguments", grid_generator_arguments);
  add_parameter("Number of refinement cycles", n_refinement_cycles);

  this->prm.enter_subsection("Error table"); // error 评估
  error_table.add_parameters(this->prm);
  this->prm.leave_subsection();
}


template <int dim> // dim模版
void
Poisson<dim>::initialize(const std::string &filename)
{
  ParameterAcceptor::initialize(filename);
}

template <int dim>
void
Poisson<dim>::parse_string(const std::string &parameters)
{
  // 方式2：将prm文件的所有的自定义参数赋到变量上面
  ParameterAcceptor::prm.parse_input_from_string(parameters);
  ParameterAcceptor::parse_all_parameters();
}



template <int dim> // dim模版
void
Poisson<dim>::make_grid()
{
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);
  triangulation.refine_global(n_refinements);
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}


template <int dim> // dim模版
void
Poisson<dim>::refine_grid()
{
  triangulation.refine_global(1);
}


template <int dim> // dim模版
void
Poisson<dim>::setup_system() // 需要增加hangding node 的处理（ppt 11）
{
  if (!fe)
    {
      fe = std::make_unique<FE_Q<dim>>(fe_degree); // fe 智能指针
      forcing_term.initialize(dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z",
                              forcing_term_expression,
                              constants);

      dirichlet_boundary_condition.initialize( //边界条件初始化赋值
        dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z",
        dirichlet_boundary_conditions_expression,
        constants);

      neumann_boundary_condition.initialize(
        dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z",
        neumann_boundary_conditions_expression,
        constants);
    }

  dof_handler.distribute_dofs(*fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  /**
   * 存储所有约束, 用于处理挂起的节点
   */
  constraints.clear();
  // create hanging node constrainints
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  for (const auto &id : dirichlet_ids)
    // will also use for boundary values from now on
    VectorTools::interpolate_boundary_values(dof_handler,
                                             id,
                                             dirichlet_boundary_condition,
                                             constraints);
  // sort, rearrange , optimise constraints
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  // Need different SparsityPattern creator
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  false); // ture的话可以移除限制
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim> // dim模版
void
Poisson<dim>::assemble_system()
{
  QGauss<dim>     quadrature_formula(fe->degree + 1);
  QGauss<dim - 1> face_quadrature_formula(
    fe->degree + 1); // 考虑边界上的积分，因此还要提前布局面积分点

  FEValues<dim> fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values |
                            update_quadrature_points);

  FEFaceValues<dim> fe_face_values(*fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_JxW_values);

  const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) +=
              (fe_values.shape_value(i, q_index) * // phi_i(x_q)
               forcing_term.value(fe_values.quadrature_point(
                 q_index)) * // f(x_q) 更新了外力项，在多阶fe节点上需要做插值
               fe_values.JxW(q_index)); // dx
        }
      // 边界上面的积分处理，类似上述的办法累加
      if (cell->at_boundary())
        //  for(const auto face: cell->face_indices())
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          if (neumann_ids.find(cell->face(f)->boundary_id()) !=
              neumann_ids.end())
            {
              fe_face_values.reinit(cell, f);
              for (const unsigned int q_index :
                   fe_face_values.quadrature_point_indices())
                for (const unsigned int i : fe_face_values.dof_indices())
                  cell_rhs(i) += fe_face_values.shape_value(i, q_index) *
                                 neumann_boundary_condition.value(
                                   fe_face_values.quadrature_point(q_index)) *
                                 fe_face_values.JxW(q_index);
            }


      cell->get_dof_indices(local_dof_indices);
      // 数据组装工作都在distribute_local_to_global完成
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  // 上面distribute_local_to_global帮你自动完成后续操作
  // std::map<types::global_dof_index, double> boundary_values;
  // VectorTools::interpolate_boundary_values(dof_handler,
  //                                          0,
  //                                          boundary_condition, //
  //                                          自定义边界条件 boundary_values);
  // MatrixTools::apply_boundary_values(boundary_values,
  //                                    system_matrix,
  //                                    solution,
  //                                    system_rhs);
}


template <int dim> // dim模版
void
Poisson<dim>::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  constraints.distribute(solution); // 最后接完需要把限制条件加回去
}


template <int dim> // dim模版
void
Poisson<dim>::output_results(const unsigned cycle) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  std::string fname =
    output_filename + "_" + std::to_string(cycle) + ".vtu"; // 定义输出的文件名
  std::ofstream output(fname);
  data_out.write_vtu(output); //输出vtu格式
}


template <int dim> // dim模版
void
Poisson<dim>::run()
{
  make_grid();
  for (unsigned int cycle = 0; cycle < n_refinement_cycles; ++cycle)
    {
      setup_system();
      assemble_system();
      solve();
      error_table.error_from_exact(
        dof_handler,
        solution,
        dirichlet_boundary_condition); // 计算相对累计误差
      output_results(cycle);           // 输出每一轮的vtu结果
      if (cycle < n_refinement_cycles - 1)
        refine_grid(); //然后根据相对误差重新细化网格
    }
  error_table.output_table(std::cout); // 最后来一个误差总表
}
