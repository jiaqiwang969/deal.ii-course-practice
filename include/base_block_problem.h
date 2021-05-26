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
#ifndef base_block_problem_include_file
#define base_block_problem_include_file

#include "base_problem.h"

// 向前声明测试者类
template <typename Integral>
class BaseBlockProblemTester;

using namespace dealii;

/**
 * 构建一个BaseBlockProblem。
 */
template <int dim>
class BaseBlockProblem : public BaseProblem<dim>
{
public:
  /**
   * Constructor. Store component names and component masks.
   */
  BaseBlockProblem(const std::vector<std::string> component_names = {{"u"}},
                   const std::string &            problem_name    = "");

  /**
   * Virtual destructor.
   */
  virtual ~BaseBlockProblem() = default;

  /**
   * 默认的CopyData对象，在WorkStream类中使用。
   */
  using CopyData = typename BaseProblem<dim>::CopyData;

  /**
   * 默认的ScratchData对象，在工作流类中使用。
   */
  using ScratchData = typename BaseProblem<dim>::ScratchData;

protected:
  /**
   * 在`cell`上组装本地系统矩阵，对FEValues和其他昂贵的抓取对象使用`scratch`，并将结果存储在`copy`对象中。
   * 关于如何使用这个函数的解释，请参见WorkStream的文档。使用这个函数。
   *
   * @param cell Cell on which we assemble the local matrix and rhs.
   * @param scratch Scratch object.
   * @param copy Copy object.
   */
  virtual void
  assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData &                                         scratch,
    CopyData &                                            copy) override;

  /**
   * 将assemble_system_on_cell()组装好的数据分配给 全局矩阵和rhs。
   *
   * @param copy 在系统矩阵和rhs上分配的本地数据。
   */
  virtual void
  copy_one_cell(const CopyData &copy) override;

  virtual void
  setup_system() override;

  virtual void
  assemble_system() override;

  virtual void
  solve() override;

  const std::vector<std::string> component_names;

  /**
   * Dofs per block
   */
  std::vector<types::global_dof_index> dofs_per_block;

  /**
   * 该MPI进程拥有的所有自由度。
   */
  std::vector<IndexSet> locally_owned_dofs;

  /**
   * 输出和误差估计所需的所有自由度。
   */
  std::vector<IndexSet> locally_relevant_dofs;

  /**
   * 系统矩阵。
   */
  LA::MPI::BlockSparseMatrix system_block_matrix;

  /**
   * 用于输出和误差估计的解决方案向量的只读副本。
   */
  LA::MPI::BlockVector locally_relevant_block_solution;

  /**
   * Solution vector.
   */
  LA::MPI::BlockVector block_solution;

  /**
   * 系统的右手边。读写向量，只包含本地拥有的dofs。
   */
  LA::MPI::BlockVector system_block_rhs;


  /**
   * 测试员类的名称。
   */
  template <typename Integral>
  friend class BaseBlockProblemTester;
};


#endif