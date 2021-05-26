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
#ifndef stokes_include_file
#define stokes_include_file

#include "base_block_problem.h"

// Forward declare the tester class
template <typename Integral>
class StokesTester;

using namespace dealii;

/**
 * 在GridGenerator命名空间中的函数可以生成的所有几何图形上，用Dirichlet或Neumann边界条件解决Poisson问题。
 */
template <int dim>
class Stokes : public BaseBlockProblem<dim>
{
public:
  /**
   * 构造函数。初始化所有参数，包括基类，并确保该类可以运行。
   */
  Stokes();

  /**
   * 销毁斯托克斯对象
   */
  ~Stokes() = default;

  using CopyData    = typename BaseBlockProblem<dim>::CopyData;
  using ScratchData = typename BaseBlockProblem<dim>::ScratchData;

protected:
  /**
   * 显示地组装斯托克斯问题。
   */
  virtual void
  assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData &                                         scratch,
    CopyData &                                            copy) override;

  virtual void
  solve() override;

  virtual void
  estimate() override;

  /**
   * 在装配程序中使用的提取器。
   */
  FEValuesExtractors::Vector velocity;

  /**
   * 在装配程序中使用的提取器。
   */
  FEValuesExtractors::Scalar pressure;

  template <typename Integral>
  friend class StokesTester;
};

#endif