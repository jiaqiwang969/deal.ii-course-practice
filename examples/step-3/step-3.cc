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
 */


// @sect3{Many new include files}

//  这些包含文件已经为你所知。它们声明了处理三角形和自由度枚举的类。
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>
// 而这是声明创建网格的函数的文件。
#include <deal.II/grid/grid_generator.h>

//  接下来的三个文件包含了所有单元的循环和从单元对象中获取信息所需的类。前两个文件以前被用来从单元中获取几何信息；最后一个是新的，它提供了一个单元中的自由度信息。
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

//  该文件包含拉格朗日插值有限元的描述。
#include <deal.II/fe/fe_q.h>

//  而这个文件是创建稀疏矩阵的稀疏模式所需要的，如前面的例子中所示。
#include <deal.II/dofs/dof_tools.h>

//  接下来的两个文件是在每个单元上使用正交法组装矩阵所需要的。下面将解释其中声明的类。
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_values.h>

//  下面三个包括我们在处理边界值时需要的文件。
#include <deal.II/base/function.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

//  我们现在几乎到了终点。第二组到最后一组include文件是用于线性代数的，我们用它来解决拉普拉斯方程的有限元离散化产生的方程组。我们将使用向量和全矩阵在每个单元中组装方程组，并将结果转移到稀疏矩阵中。然后我们将使用共轭梯度求解器来解决这个问题，为此我们需要一个预处理程序（在这个程序中，我们使用身份预处理程序，它没有任何作用，但我们还是需要包含这个文件）。
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

// 最后，这是用来输出到文件和控制台的。
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

// ...这是为了将deal.II命名空间导入到全局范围。
using namespace dealii;

// @sect3{The <code>Step3</code> class}


//  在这个程序中，我们没有采用以前例子中的程序化编程，而是将所有东西都封装在一个类中。该类由各自执行有限元程序某些方面的函数组成，一个控制先做什么和后做什么的`main`函数，以及一个成员变量的列表。

//  该类的公共部分相当短：它有一个构造函数和一个从外部调用的函数`run`，其作用类似于`main`函数：它协调该类的哪些操作应以何种顺序运行。类中的其他东西，也就是所有真正做事情的函数，都在类的私有部分。
class Step3
{
public:
  Step3();

  void
  run();

  //  然后是成员函数，它们主要做它们名字所暗示的事情，在介绍中已经讨论过了。因为它们不需要从外部调用，所以它们是这个类的私有函数。

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results() const;

  // 最后我们还有一些成员变量。有一些变量描述了三角形和自由度的全局编号（我们将在这个类的构造函数中指定有限元的确切多项式程度）...
  Triangulation<2> triangulation;
  FE_Q<2>          fe;
  DoFHandler<2>    dof_handler;

  //  ...拉普拉斯方程离散化产生的系统矩阵的稀疏模式和数值的变量...
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  //  ...和变量，这些变量将持有右手边和解决方案的向量。
  Vector<double> solution;
  Vector<double> system_rhs;
};

// @sect4{Step3::Step3}

//  这里有一个构造函数。它除了首先指定我们需要双线性元素（由有限元对象的参数表示，它表示多项式的程度），并将dof_handler变量与我们使用的三角形相关联之外，没有做更多的工作。(注意，目前三角结构并没有设置网格，但是DoFHandler并不关心：它只想知道它将与哪个三角结构相关联，只有当你使用distribution_dofs()函数试图在网格上分布自由度时，它才开始关心实际的网格。)
//  Step3类的所有其他成员变量都有一个默认的构造函数，它可以完成我们想要的一切。
Step3::Step3()
  : fe(1)
  , dof_handler(triangulation)
{}


// @sect4{Step3::make_grid}

//  现在，我们要做的第一件事是生成我们想在其上进行计算的三角形，并对每个顶点进行自由度编号。我们之前在步骤1和步骤2中分别看到了这两个步骤。

//  这个函数做的是第一部分，创建网格。
//  我们创建网格并对所有单元格进行五次细化。由于初始网格（即正方形$[-1,1]\times[-1,1]$）只由一个单元组成，最终的网格有32乘以32个单元，总共1024个。

//  不确定1024是否是正确的数字？我们可以通过使用三角形上的<code>n_active_cells()</code>函数输出单元格的数量来检查。
void
Step3::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}

//  @note 我们调用Triangulation::n_active_cells()函数，而不是Triangulation::n_cells()。这里，<i>active</i>指的是没有被进一步细化的单元。我们强调 "活跃 "这个形容词，因为还有更多的单元，即最细单元的父单元，它们的父单元等等，直到构成初始网格的一个单元。当然，在下一个更粗的层次上，单元格的数量是最细层次上单元格的四分之一，即256，然后是64、16、4和1。如果你在上面的代码中调用<code>triangulation.n_cells()</code>，你会因此得到一个1365的值。另一方面，单元格的数量（相对于活动单元格的数量）通常没有什么意义，所以没有很好的理由去打印它。


// @sect4{Step3::setup_system}

//  接下来我们列举所有的自由度，并设置矩阵和向量对象来保存系统数据。枚举是通过使用DoFHandler::distribution_dofs()完成的，正如我们在步骤2的例子中看到的那样。由于我们使用了FE_Q类，并且在构造函数中设置了多项式的度数为1，即双线性元素，这就将一个自由度与每个顶点联系起来。当我们在生成输出时，让我们也看看有多少自由度被生成。
void
Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  //  每个顶点应该有一个DoF。因为我们有一个32乘以32的网格，所以DoF的数量应该是33乘以33，即1089。

  //  正如我们在前面的例子中所看到的，我们通过首先创建一个临时结构，标记那些可能为非零的条目，然后将数据复制到SparsityPattern对象中，然后可以被系统矩阵使用，来设置一个稀疏模式。
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  //  请注意，SparsityPattern对象并不持有矩阵的值，它只存储条目所在的位置。条目本身被存储在SparseMatrix类型的对象中，我们的变量system_matrix就是其中之一。区分稀疏模式和矩阵是为了让几个矩阵使用相同的稀疏模式。这在这里似乎并不重要，但是当你考虑到矩阵的大小，以及建立稀疏模式可能需要一些时间时，如果你必须在程序中存储几个矩阵，这在大规模问题中就变得很重要。
  system_matrix.reinit(sparsity_pattern);

  //  在这个函数中，最后要做的是将右侧向量和解向量的大小设置为正确的值。
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

// @sect4{Step3::assemble_system}


//  下一步是计算形成线性系统的矩阵和右手边的条目，我们从中计算出解决方案。这是每个有限元程序的核心功能，我们在介绍中已经讨论了主要步骤。

//  组装矩阵和向量的一般方法是在所有单元上循环，并在每个单元上通过正交计算该单元对全局矩阵和右侧的贡献。现在要认识到的一点是，我们需要实心单元上正交点位置的形状函数值。然而，有限元形状函数和正交点都只定义在参考单元上。因此，它们对我们帮助不大，事实上，我们几乎不会直接从这些对象中查询有关有限元形状函数或正交点的信息。

//  相反，我们需要的是一种将这些数据从参考单元映射到实际单元的方法。可以做到这一点的类是由Mapping类派生出来的，尽管人们常常不必直接与它们打交道：库中的许多函数可以将映射对象作为参数，但当它被省略时，它们只是求助于标准的双线性Q1映射。我们将走这条路，暂时不去管它（我们在第10步、第11步和第12步再来讨论这个问题）。

//  所以我们现在有三个类的集合来处理：有限元、正交和映射对象。这就太多了，所以有一种类型的类可以协调这三者之间的信息交流：FEValues类。如果给这三个对象各一个实例（或两个，以及一个隐式线性映射），它就能为你提供实心单元上正交点的形状函数的值和梯度的信息。

//  利用所有这些，我们将把这个问题的线性系统组装在以下函数中。
void
Step3::assemble_system()
{
  //  好吧，让我们开始吧：我们需要一个正交公式来评估每个单元的积分。让我们采用一个高斯公式，每个方向有两个正交点，即总共有四个点，因为我们是在二维。这个正交公式可以精确地积分三度以下的多项式（在一维）。很容易检查出，这对目前的问题来说是足够的。
  QGauss<2> quadrature_formula(fe.degree + 1);
  //  然后我们初始化我们在上面简单谈过的对象。它需要被告知我们要使用哪个有限元，以及正交点和它们的权重（由一个正交对象共同描述）。如前所述，我们使用隐含的Q1映射，而不是自己明确指定一个。最后，我们必须告诉它我们希望它在每个单元上计算什么：我们需要正交点的形状函数值（对于右手边的$(\varphi_i,f)$），它们的梯度（对于矩阵条目$(\nabla
  //  \varphi_i, \nabla
  //  \varphi_j)$），还有正交点的权重和从参考单元到实际单元的雅各布变换的行列式。

  //  我们实际需要的信息清单是作为FEValues构造函数的第三个参数，以标志集合的形式给出的。由于这些值必须重新计算，或更新，每次我们去一个新的单元，所有这些标志都以<code>update_</code>的前缀开始，然后指出我们想要更新的实际内容。如果我们想要计算形状函数的值，那么给出的标志是#update_values；对于梯度，它是#update_gradients。雅各布的行列式和正交权重总是一起使用的，所以只有乘积（雅各布乘以权重，或简称<code>JxW</code>）被计算；因为我们需要它们，所以我们也必须列出#update_JxW_values。
  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);
  //  这种方法的优点是，我们可以指定在每个单元上实际需要什么样的信息。很容易理解的是，这种方法可以大大加快有限元计算的速度，相比之下，所有的东西，包括二阶导数、单元的法向量等都在每个单元上计算，不管是否需要它们。

  //  @注意 <code>update_values | update_gradients | update_JxW_values</code> 语法对于那些不习惯用C语言编程多年的位操作的人来说不是很明显。首先，<code>operator|</code>是<i>bitwise or operator</i>，也就是说，它接受两个整数参数，这些参数被解释为位模式，并返回一个整数，其中每个位都被设置，因为在两个参数中至少有一个的对应位被设置。例如，考虑操作<code>9|10</code>。在二进制中，<code>9=0b1001</code>（其中前缀<code>0b</code>表示该数字将被解释为二进制数字）和<code>10=0b1010</code>。通过每个比特，看它是否在其中一个参数中被设置，我们得出<code>0b1001|0b1010=0b1011</code>，或者，用十进制符号表示，<code>9|10=11</code>。你需要知道的第二个信息是，各种<code>update_*</code>标志都是整数，<i>正好有一个比特被设置</i>。例如，假设<code>update_values=0b00001=1</code>，<code>update_gradients=0b00010=2</code>，<code>update_JxW_values=0b10000=16</code>。那么<code>update_values | update_gradients | update_JxW_values = 0b10011 = 19</code>。换句话说，我们得到一个数字，<i>编码一个二进制掩码，代表你想要发生的所有操作</i>，其中每个操作正好对应于整数中的一个位，如果等于1，意味着每个单元格上的特定片断应该被更新，如果是0，意味着我们不需要计算它。换句话说，即使<code>operator|</code>是<i>bitwise OR操作</i>，它真正代表的是<i>我想要这个和那个和另一个</i>。这样的二进制掩码在C语言编程中很常见，但在C++这样的高级语言中也许不是这样，但却能很好地满足当前的目的。

  //  为了在下文中进一步使用，我们为一个将被频繁使用的值定义了一个快捷方式。也就是每个单元的自由度数量的缩写（因为我们是在二维，自由度只与顶点相关，所以这个数字是四，但是我们更希望在写这个变量的定义时，不妨碍我们以后选择不同的有限元，每个单元有不同的自由度数量，或者在不同的空间维度工作）。

  //  一般来说，使用符号名称而不是硬编码这些数字是个好主意，即使你知道它们，因为例如，你可能想在某个时候改变有限元。改变元素就必须在不同的函数中进行，而且很容易忘记在程序的另一部分做相应的改变。最好不要依赖自己的计算，而是向正确的对象索取信息。在这里，我们要求有限元告诉我们每个单元的自由度数，不管我们在程序的其他地方选择的空间维度或多项式程度如何，我们都会得到正确的数字。

  //  这里定义的快捷方式主要是为了讨论基本概念，而不是因为它可以节省大量的打字量，那么下面的循环就会变得更容易阅读。在大型程序中，你会在很多地方看到这样的快捷方式，`dofs_per_cell`就是一个或多或少是这类对象的常规名称。
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  // Now, we said that we wanted to assemble the global matrix and vector
  // cell-by-cell. We could write the results directly into the global matrix,
  // but this is not very efficient since access to the elements of a sparse
  // matrix is slow. Rather, we first compute the contribution of each cell in
  // a small matrix with the degrees of freedom on the present cell, and only
  // transfer them to the global matrix when the computations are finished for
  // this cell. We do the same for the right hand side vector. So let's first
  // allocate these objects (these being local objects, all degrees of freedom
  // are coupling with all others, and we should use a full matrix object
  // rather than a sparse one for the local operations; everything will be
  // transferred to a global sparse matrix later on):
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // When assembling the contributions of each cell, we do this with the local
  // numbering of the degrees of freedom (i.e. the number running from zero
  // through dofs_per_cell-1). However, when we transfer the result into the
  // global matrix, we have to know the global numbers of the degrees of
  // freedom. When we query them, we need a scratch (temporary) array for
  // these numbers (see the discussion at the end of the introduction for
  // the type, types::global_dof_index, used here):
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Now for the loop over all cells. We have seen before how this works for a
  // triangulation. A DoFHandler has cell iterators that are exactly analogous
  // to those of a Triangulation, but with extra information about the degrees
  // of freedom for the finite element you're using. Looping over the active
  // cells of a degree-of-freedom handler works the same as for a triangulation.
  //
  // Note that we declare the type of the cell as `const auto &` instead of
  // `auto` this time around. In step 1, we were modifying the cells of the
  // triangulation by flagging them with refinement indicators. Here we're only
  // examining the cells without modifying them, so it's good practice to
  // declare `cell` as `const` in order to enforce this invariant.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // We are now sitting on one cell, and we would like the values and
      // gradients of the shape functions be computed, as well as the
      // determinants of the Jacobian matrices of the mapping between
      // reference cell and true cell, at the quadrature points. Since all
      // these values depend on the geometry of the cell, we have to have the
      // FEValues object re-compute them on each cell:
      fe_values.reinit(cell);

      // Next, reset the local cell's contributions to global matrix and
      // global right hand side to zero, before we fill them:
      cell_matrix = 0;
      cell_rhs    = 0;

      // Now it is time to start integration over the cell, which we
      // do by looping over all quadrature points, which we will
      // number by q_index.
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          // First assemble the matrix: For the Laplace problem, the
          // matrix on each cell is the integral over the gradients of
          // shape function i and j. Since we do not integrate, but
          // rather use quadrature, this is the sum over all
          // quadrature points of the integrands times the determinant
          // of the Jacobian matrix at the quadrature point times the
          // weight of this quadrature point. You can get the gradient
          // of shape function $i$ at quadrature point with number q_index by
          // using <code>fe_values.shape_grad(i,q_index)</code>; this
          // gradient is a 2-dimensional vector (in fact it is of type
          // Tensor@<1,dim@>, with here dim=2) and the product of two
          // such vectors is the scalar product, i.e. the product of
          // the two shape_grad function calls is the dot
          // product. This is in turn multiplied by the Jacobian
          // determinant and the quadrature point weight (that one
          // gets together by the call to FEValues::JxW() ). Finally,
          // this is repeated for all shape functions $i$ and $j$:
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx

          // We then do the same thing for the right hand side. Here,
          // the integral is over the shape function i times the right
          // hand side function, which we choose to be the function
          // with constant value one (more interesting examples will
          // be considered in the following programs).
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            1 *                                 // f(x_q)
                            fe_values.JxW(q_index));            // dx
        }
      // Now that we have the contribution of this cell, we have to transfer
      // it to the global matrix and right hand side. To this end, we first
      // have to find out which global numbers the degrees of freedom on this
      // cell have. Let's simply ask the cell for that information:
      cell->get_dof_indices(local_dof_indices);

      // Then again loop over all shape functions i and j and transfer the
      // local elements to the global matrix. The global numbers can be
      // obtained using local_dof_indices[i]:
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

      // And again, we do the same thing for the right hand side vector.
      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  // Now almost everything is set up for the solution of the discrete
  // system. However, we have not yet taken care of boundary values (in fact,
  // Laplace's equation without Dirichlet boundary values is not even uniquely
  // solvable, since you can add an arbitrary constant to the discrete
  // solution). We therefore have to do something about the situation.
  //
  // For this, we first obtain a list of the degrees of freedom on the
  // boundary and the value the shape function shall have there. For
  // simplicity, we only interpolate the boundary value function, rather than
  // projecting it onto the boundary. There is a function in the library which
  // does exactly this: VectorTools::interpolate_boundary_values(). Its
  // parameters are (omitting parameters for which default values exist and
  // that we don't care about): the DoFHandler object to get the global
  // numbers of the degrees of freedom on the boundary; the component of the
  // boundary where the boundary values shall be interpolated; the boundary
  // value function itself; and the output object.
  //
  // The component of the boundary is meant as follows: in many cases, you may
  // want to impose certain boundary values only on parts of the boundary. For
  // example, you may have inflow and outflow boundaries in fluid dynamics, or
  // clamped and free parts of bodies in deformation computations of
  // bodies. Then you will want to denote these different parts of the
  // boundary by indicators, and tell the interpolate_boundary_values
  // function to only compute the boundary values on a certain part of the
  // boundary (e.g. the clamped part, or the inflow boundary). By default,
  // all boundaries have a 0 boundary indicator, unless otherwise specified. If
  // sections of the boundary have different boundary conditions, you have to
  // number those parts with different boundary indicators. The function call
  // below will then only determine boundary values for those parts of the
  // boundary for which the boundary indicator is in fact the zero specified as
  // the second argument.
  //
  // The function describing the boundary values is an object of type Function
  // or of a derived class. One of the derived classes is
  // Functions::ZeroFunction, which describes (not unexpectedly) a function
  // which is zero everywhere. We create such an object in-place and pass it to
  // the VectorTools::interpolate_boundary_values() function.
  //
  // Finally, the output object is a list of pairs of global degree of freedom
  // numbers (i.e. the number of the degrees of freedom on the boundary) and
  // their boundary values (which are zero here for all entries). This mapping
  // of DoF numbers to boundary values is done by the <code>std::map</code>
  // class.
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<2>(),
                                           boundary_values);
  // Now that we got the list of boundary DoFs and their respective boundary
  // values, let's use them to modify the system of equations
  // accordingly. This is done by the following function call:
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


// @sect4{Step3::solve}

// The following function simply solves the discretized equation. As the
// system is quite a large one for direct solvers such as Gauss elimination or
// LU decomposition, we use a Conjugate Gradient algorithm. You should
// remember that the number of variables here (only 1089) is a very small
// number for finite element computations, where 100.000 is a more usual
// number.  For this number of variables, direct methods are no longer usable
// and you are forced to use methods like CG.
void
Step3::solve()
{
  // First, we need to have an object that knows how to tell the CG algorithm
  // when to stop. This is done by using a SolverControl object, and as
  // stopping criterion we say: stop after a maximum of 1000 iterations (which
  // is far more than is needed for 1089 variables; see the results section to
  // find out how many were really used), and stop if the norm of the residual
  // is below $10^{-12}$. In practice, the latter criterion will be the one
  // which stops the iteration:
  SolverControl solver_control(1000, 1e-12);
  // Then we need the solver itself. The template parameter to the SolverCG
  // class is the type of the vectors, but the empty angle brackets indicate
  // that we simply take the default argument (which is
  // <code>Vector@<double@></code>):
  SolverCG<Vector<double>> solver(solver_control);

  // Now solve the system of equations. The CG solver takes a preconditioner
  // as its fourth argument. We don't feel ready to delve into this yet, so we
  // tell it to use the identity operation as preconditioner:
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  // Now that the solver has done its job, the solution variable contains the
  // nodal values of the solution function.
}


// @sect4{Step3::output_results}

// The last part of a typical finite element program is to output the results
// and maybe do some postprocessing (for example compute the maximal stress
// values at the boundary, or the average flux across the outflow, etc). We
// have no such postprocessing here, but we would like to write the solution
// to a file.
void
Step3::output_results() const
{
  // To write the output to a file, we need an object which knows about output
  // formats and the like. This is the DataOut class, and we need an object of
  // that type:
  DataOut<2> data_out;
  // Now we have to tell it where to take the values from which it shall
  // write. We tell it which DoFHandler object to use, and the solution vector
  // (and the name by which the solution variable shall appear in the output
  // file). If we had more than one vector which we would like to look at in
  // the output (for example right hand sides, errors per cell, etc) we would
  // add them as well:
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  // After the DataOut object knows which data it is to work on, we have to
  // tell it to process them into something the back ends can handle. The
  // reason is that we have separated the frontend (which knows about how to
  // treat DoFHandler objects and data vectors) from the back end (which knows
  // many different output formats) and use an intermediate data format to
  // transfer data from the front- to the backend. The data is transformed
  // into this intermediate format by the following function:
  data_out.build_patches();

  // Now we have everything in place for the actual output. Just open a file
  // and write the data into it, using VTK format (there are many other
  // functions in the DataOut class we are using here that can write the
  // data in postscript, AVS, GMV, Gnuplot, or some other file
  // formats):
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
}


// @sect4{Step3::run}

// Finally, the last function of this class is the main function which calls
// all the other functions of the <code>Step3</code> class. The order in which
// this is done resembles the order in which most finite element programs
// work. Since the names are mostly self-explanatory, there is not much to
// comment about:
void
Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}


// @sect3{The <code>main</code> function}

// This is the main function of the program. Since the concept of a
// main function is mostly a remnant from the pre-object oriented era
// before C++ programming, it often does not do much more than
// creating an object of the top-level class and calling its principle
// function.
//
// Finally, the first line of the function is used to enable output of
// some diagnostics that deal.II can generate.  The @p deallog
// variable (which stands for deal-log, not de-allog) represents a
// stream to which some parts of the library write output. For
// example, iterative solvers will generate diagnostics (starting
// residual, number of solver steps, final residual) as can be seen
// when running this tutorial program.
//
// The output of @p deallog can be written to the console, to a file,
// or both. Both are disabled by default since over the years we have
// learned that a program should only generate output when a user
// explicitly asks for it. But this can be changed, and to explain how
// this can be done, we need to explain how @p deallog works: When
// individual parts of the library want to log output, they open a
// "context" or "section" into which this output will be placed. At
// the end of the part that wants to write output, one exits this
// section again. Since a function may call another one from within
// the scope where this output section is open, output may in fact be
// nested hierarchically into these sections. The LogStream class of
// which @p deallog is a variable calls each of these sections a
// "prefix" because all output is printed with this prefix at the left
// end of the line, with prefixes separated by colons. There is always
// a default prefix called "DEAL" (a hint at deal.II's history as the
// successor of a previous library called "DEAL" and from which the
// LogStream class is one of the few pieces of code that were taken
// into deal.II).
//
// By default, @p logstream only outputs lines with zero prefixes --
// i.e., all output is disabled because the default "DEAL" prefix is
// always there. But one can set a different maximal number of
// prefixes for lines that should be output to something larger, and
// indeed here we set it to two by calling
// LogStream::depth_console(). This means that for all screen output,
// a context that has pushed one additional prefix beyond the default
// "DEAL" is allowed to print its output to the screen ("console"),
// whereas all further nested sections that would have three or more
// prefixes active would write to @p deallog, but @p deallog does not
// forward this output to the screen. Thus, running this example (or
// looking at the "Results" section), you will see the solver
// statistics prefixed with "DEAL:CG", which is two prefixes. This is
// sufficient for the context of the current program, but you will see
// examples later on (e.g., in step-22) where solvers are nested more
// deeply and where you may get useful information by setting the
// depth even higher.
int
main()
{
  deallog.depth_console(2);

  Step3 laplace_problem;
  laplace_problem.run();

  return 0;
}
