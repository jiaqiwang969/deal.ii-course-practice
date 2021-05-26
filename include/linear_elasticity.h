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
#ifndef linear_elasticity_include_file
#define linear_elasticity_include_file

#include "base_problem.h"

// Forward declare the tester class
template <typename Integral>
class LinearElasticityTester;

using namespace dealii;

/**
 * 在所有可由GridGenerator中的函数生成的几何图形上，用Dirichlet或Neumann边界条件解决泊松问题。
 * 名称空间中的函数生成的所有几何图形。
 */
template <int dim>
class LinearElasticity : public BaseProblem<dim>
{
public:
  /**
   *  构造函数。初始化所有参数，包括基类，并确保该类可以运行。
   */
  LinearElasticity();

  /**
   * 销毁线性弹性对象
   */
  ~LinearElasticity() = default;

  using CopyData    = typename BaseProblem<dim>::CopyData;
  using ScratchData = typename BaseProblem<dim>::ScratchData;

protected:
  /**
   * 显式组装线性弹性问题。
   */
  virtual void
  assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData &                                         scratch,
    CopyData &                                            copy);

  /**
   * 在装配程序中使用的提取器。
   */
  FEValuesExtractors::Vector velocity;

  // 线性弹性的物理参数
  double mu     = 1;
  double lambda = 1;

  template <typename Integral>
  friend class LinearElasticityTester;
};

#endif
