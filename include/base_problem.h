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
#ifndef base_problem_include_file
#define base_problem_include_file

#include <deal.II/base/function.h> // 提供了一些零函数、常函数
#include <deal.II/base/function_parser.h> // 函数转化为代码，FunctionParser
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parameter_acceptor.h> // 参数接收，以及输送到prm
#include <deal.II/base/parsed_convergence_table.h> // 结果分析整理成table
#include <deal.II/base/quadrature_lib.h>           // 不同的积分点策略
#include <deal.II/base/thread_management.h>
#include <deal.II/base/timer.h> // 计时
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/grid_refinement.h> // 分布式网格加密策略
#include <deal.II/distributed/tria.h>            // 分布式网格划分策略

#include <deal.II/dofs/dof_handler.h> // 自由度分配管理策略
#include <deal.II/dofs/dof_tools.h> //自由度相关的工具，比如make_sparsity_pattern

#include <deal.II/fe/fe_q.h> // this->degree
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h> // 有限元配置，积分点，mapping的一次大装配
#include <deal.II/fe/fe_values_extractors.h> // 允许你将单一的形状函数解释为张量、标量等类型的对象
#include <deal.II/fe/mapping_q_generic.h> // 导数矩阵，插值矩阵等，利用映射关系

#include <deal.II/grid/grid_generator.h> // 基本的形状网格，cube等
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h> // 网格划分

#include <deal.II/lac/affine_constraints.h> // 此处用于处理hanging nodes节点的不连续问题
#include <deal.II/lac/dynamic_sparsity_pattern.h> // 动态稀疏矩阵存储，避免内存过大
#include <deal.II/lac/full_matrix.h> // 矩阵的一些相关操作，例如求转置等
#include <deal.II/lac/generic_linear_algebra.h> // 分布式下的petsc或 Trilinos预解方程
#include <deal.II/lac/precondition.h> // 预解方程
#include <deal.II/lac/solver_cg.h>    // CG算法求解方程
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h> // 稀疏矩阵相关算法
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/vector.h> // 向量相关

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/copy_data.h> // 拷贝相关的数据操作
#include <deal.II/meshworker/scratch_data.h>
#include <deal.II/meshworker/scratch_data.h> // 消耗cpu的计算操作

#include <deal.II/numerics/data_out.h> // 数据后处理相关操作
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h> // 处理矩阵的相关类库
#include <deal.II/numerics/vector_tools.h> // 处理向量的相关类库

#include <boost/signals2.hpp> // 插眼

#include <fstream>  // 文件流，跟文件处理相关的操作
#include <iostream> // 字符串相关操作


/**
 * 选择分布式矩阵计算外接库：PETSC或者TRILINOS
 */
#define FORCE_USE_OF_TRILINOS

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA


// 声明测试类
template <typename Integral>
class BaseProblemTester;

using namespace dealii;

/**
 * Solve the Poisson problem, with Dirichlet or Neumann boundary conditions, on
 * all geometries that can be generated by the functions in the GridGenerator
 * namespace.
 */
template <int dim>
class BaseProblem : public ParameterAcceptor // 接受参数的传递
{
public:
  /**
   * 构造函数。初始化所有参数，并确保该类可以运行。
   */
  BaseProblem(const unsigned int &n_components = 1,
              const std::string & problem_name = ""); // n_components 默认为1

  /**
   * Virtual destructor.
   */
  virtual ~BaseProblem() = default; // 析构函数

  /**
   * 输出一些琐碎的信息，如dofs的数量、单元格、线程等。
   */
  void
  print_system_info();


  /**
   * 问题的主要切入点。
   */
  void
  run();


  /**
   * 用给定的文件初始化内部参数。
   *
   * @param filename 指的是参数的文件名.
   */
  void
  initialize(const std::string &filename);


  /**
   * 解析一个字符串，好像它是一个参数文件，并相应地设置参数。主要用在测试中，方便测试过程。
   *
   * @param par 包含一些参数的字符串。
   */
  void
  parse_string(const std::string &par); // 为了直接从test中嵌入代码,进行相关测试

  /**
   * 默认的CopyData对象，在WorkStream类中使用。
   */
  using CopyData = MeshWorker::CopyData<1, 1, 1>;
  /**
   * 默认的ScratchData对象，在工作流类中使用。
   */
  using ScratchData = MeshWorker::ScratchData<dim>;

protected:
  /**
   * 在`cell`上组装本地系统矩阵，对FEValues和其他昂贵的抓取对象使用`scratch`，并将结果存储在`copy`对象中。
   * 关于如何使用这个函数的解释，请看WorkStream的文档。
   *
   * @param cell 我们在其上组装本地矩阵和rhs的单元。
   * @param scratch Scratch object.
   * @param copy Copy object.
   */
  virtual void // 虚函数：为了允许用基类的指针来调用子类的这个函数
  assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData &                                         scratch,
    CopyData &                                            copy);


  /**
   * 将assemble_system_on_cell()组装好的数据分配给全局矩阵和rhs。
   *
   * @param copy 在系统矩阵和rhs上分配的本地数据。
   */
  virtual void
  copy_one_cell(const CopyData &copy);


  /**
   * 生成参数文件中指定的初始网格。
   */
  void
  make_grid();

  /**
   * 求解全局系统。
   */
  virtual void
  solve();


  /**
   *  进行后验误差估计，并将结果存储在`error_per_cell`向量中。
   */
  virtual void
  estimate();


  /**
   * 根据所选择的策略，标记一些单元格进行细化。
   */
  void
  mark();


  /**
   * 细化网格
   */
  void
  refine_grid();

  /**
   *   初始设置：分配自由度，使所有向量和矩阵的大小合适，初始化函数和指针。
   */
  virtual void
  setup_system();

  /**
   * 在setup_system()结束时调用的信号。
   */
  boost::signals2::signal<void()> setup_system_call_back;

  /**
   * 连接到这个信号来增加数据向量。
   */
  boost::signals2::signal<void(DataOut<dim> &)> add_data_vector;

  /**
   * 实际上是在单元格上循环，并组装全局系统。
   */
  virtual void
  assemble_system();


  /**
   * 以Paraview或Visit可以读取的格式输出解决方案和网格。
   *
   * @param cycle 网格加密循环次数
   */
  virtual void
  output_results(const unsigned cycle) const;


  /**
   * 组份的数目
   */
  const unsigned int n_components;

  /**
   * 全局 mpi 通信.
   */
  MPI_Comm mpi_communicator;

  /**
   * 用于多线程装配的线程数。
   */
  int number_of_threads = -1;

  /**
   * 只在零号处理器上输出。
   */
  ConditionalOStream pcout;

  /**
   * 时间信息输出。
   */
  mutable TimerOutput timer;

  /**
   * 网格划分策略。
   */
  parallel::distributed::Triangulation<dim> triangulation;


  /**
   * 有限元空间
   *
   * 这是一个独特的指针，允许通过参数文件创建。
   */
  std::unique_ptr<FiniteElement<dim>> fe;

  /**
   * 参考元素和真实元素之间的映射。
   *
   * 这是一个独特的指针，允许通过参数文件创建。
   */
  std::unique_ptr<MappingQGeneric<dim>> mapping;

  /**
   * 自由度的处理者。
   */
  DoFHandler<dim> dof_handler;

  /**
   * 悬挂的节点和基本的边界条件。
   */
  AffineConstraints<double> constraints;


  /**
   * 该MPI进程拥有的所有自由度。
   */
  IndexSet locally_owned_dofs;

  /**
   * 输出和误差估计所需的所有自由度。
   */
  IndexSet locally_relevant_dofs;

  /**
   * 系统矩阵
   */
  LA::MPI::SparseMatrix system_matrix;

  /**
   * 用于输出和误差估计的向量解的只读副本。
   */
  LA::MPI::Vector locally_relevant_solution;

  /**
   * 向量解。
   */
  LA::MPI::Vector solution;

  /**
   * 系统的右手边。读写向量，只包含本地拥有的自由度。
   */
  LA::MPI::Vector system_rhs;


  /**
   * 本地误差估计器的存储。这个向量也包含与人工单元相关的值（即它的长度为`triangulation.n_active_cells()`），
   * 但它只在本地拥有的单元上非零。estimate()方法只填入本地拥有的单元。
   */
  Vector<float> error_per_cell;

  /**
   * 在 "精确 "估计器、"凯利 "估计器和 "残差 "估计器之间选择。
   */
  std::string estimator_type = "exact";

  /**
   * 在 "全局"、"固定_分数"（也称为Dorfler标记策略）和 "固定_数字"之间选择。
   */
  std::string marking_strategy = "global";

  /**
   * 粗化和细化的分数。
   */
  std::pair<double, double> coarsening_and_refinement_factors = {0.03, 0.3};

  /**
   * 任何直径小于该函数在单元中心评估的单元将被再次加密。这在开始模拟之前要进行`n_refimement`次。
   */
  FunctionParser<dim> pre_refinement;

  /**
   * 外力项放在方程的右边。
   */
  FunctionParser<dim> forcing_term;

  /**
   * 用来计算误差的函数。
   */
  FunctionParser<dim> exact_solution;

  /**
   *非均质dirichlet边界条件。
   */
  FunctionParser<dim> dirichlet_boundary_condition;

  /**
   * 非均质neumann边界条件。
   */
  FunctionParser<dim> neumann_boundary_condition;

  /**
   * 用来生成有限元空间的字符串。应该与FETools::get_fe_by_name()兼容。
   */
  std::string fe_name = "FE_Q(1)";

  /**
   * 参考元素和实际元素之间的映射程度。
   */
  unsigned int mapping_degree = 1;

  /**
   * 仿真开始前要进行的预加密的次数。
   */
  unsigned int n_refinements = 4;

  /**
   * 要执行的求解-估计-标记-细化循环的数量。
   */
  unsigned int n_refinement_cycles = 1;

  /**
   * 解决方案的输出名称。
   */
  std::string output_filename = "linear_elasticity";

  /**
   * 在哪些边界id上，我们施加基本的边界条件。
   */
  std::set<types::boundary_id> dirichlet_ids = {0};

  /**
   * 在哪些边界id上，我们施加基本的边界条件。
   */
  std::set<types::boundary_id> neumann_ids;

  /**
   * 在各种函数对象的定义中使用的常量。
   */
  std::map<std::string, double> constants;

  /**
   * 外力项的表达式。
   */
  std::string forcing_term_expression = "1";

  /**
   * 精确解的表达式
   */
  std::string exact_solution_expression = "0";

  /**
   * Dirichlet 边界条件的表达式。
   */
  std::string dirichlet_boundary_conditions_expression = "0";

  /**
   * Neumann 边界条件的表达式。
   */
  std::string neumann_boundary_conditions_expression = "0";

  /**
   * pre-refinement 的表达式。
   */
  std::string pre_refinement_expression = "0";

  /**
   * 要调用的GridGenerator函数的名称。
   */
  std::string grid_generator_function = "hyper_cube";

  /**
   * 传递给GridGenerator函数的参数。
   */
  std::string grid_generator_arguments = "0: 1: false";

  /**
   * 一个用于输出收敛错误的表格。
   */
  ParsedConvergenceTable error_table;

  /**
   * 用于存储求解器参数的类，如最大迭代次数、绝对公差和相对公差。
   */
  ParameterAcceptorProxy<ReductionControl> solver_control;

  /**
   * 测试员类的名称。
   */
  template <typename Integral>
  friend class BaseProblemTester;
};


/**
 * 基本的流水线过程，也凝练成了类
 */
template <typename ProblemType>
int
run(int argc, char **argv)
{
  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
      std::string                      par_name = "";
      if (argc > 1)
        par_name = argv[1];

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog.depth_console(2);
      else
        deallog.depth_console(0);

      ProblemType base_problem;
      base_problem.initialize(par_name);
      base_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}


#endif
