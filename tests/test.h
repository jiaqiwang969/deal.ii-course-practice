//*/// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


/**







   @page DEALGlossary Glossary



 本词汇表解释了一些在deal.II类的文档中经常使用的术语。词汇表通常只给出了一个特定概念的微观视图；如果你对大局感到困惑，那么也值得参考 @ref index 页面上的类的全局概述。



   <dl>



 <dt class="glossary">  @anchor  GlossActive <b>Active cells</b></dt>  <dd>
如果一个单元格、面或边没有被进一步细化，即没有子代，则被定义为<i>active</i>
。一旦一个单元、面或边成为父级，它就不再活跃。除非使用多网格算法，否则活动单元是唯一携带自由度的单元。
</dd>











 <dt class="glossary">  @anchor  GlossArtificialCell <b>Artificial cells</b></dt>  <dd>  如果一个网格使用 parallel::distributed::Triangulation 类分布在多个MPI进程中，每个处理器只存储它拥有的单元，一层由其他处理器拥有的相邻单元（称为 @ref GlossGhostCell  "幽灵单元"），所有粗级单元，以及所有为保持相邻单元必须最多有一个细化级别这一不变性所需的单元。存储在每个进程上的不属于该进程且不属于幽灵单元的单元被称为 "人工单元"，对于这些单元，谓词 <code>cell-@>is_artificial()</code> 返回真。人工单元被保证存在于全局分布的网格中，但它们可能在其他处理器上被进一步细化。更多信息请参见 @ref distributed_paper  "分布式计算论文"。



 人工单元的概念对于在每个处理器上存储整个网格的三角计算没有意义，即
dealii::Triangulation 类。   </dd>







 <dt class="glossary">  @anchor  GlossBlockLA <b>Block (linear
algebra)</b></dt>。



   <dd>  将一个矩阵或向量视为单个块的集合通常很方便。例如，在 step-20
（和其他教程程序）中，我们要考虑全局线性系统 $Ax=b$ 的形式@f{eqnarray*}
  \left(\begin{array}{cc}
    M & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),
 @f}，其中 $U,P$ 分别是速度和压力自由度的值， $M$ 是速度空间上的质量矩阵， $B$
对应于负发散算子， $B^T$ 是其转位，对应于负梯度。



 使用这种分解为块的方法，人们可以定义基于方程组中存在的单个算子（例如，在
step-20
的情况下，Schur补码）的预处理程序，而不是整个矩阵。实质上，块被用来反映线性代数中PDE系统的结构，特别是允许对具有多个解决方案组件的问题进行模块化求解。另一方面，矩阵和右手边的向量也可以作为一个单元来处理，这在线性系统的装配过程中是很方便的，因为人们可能不想对各个组件进行区分，或者对于不关心块结构的外Krylov空间求解器（例如，如果只有预处理程序需要块结构）。



 将矩阵和向量分割成块是由BlockSparseMatrix、BlockVector和相关类支持的。参见 @ref LAC 模块中对各种线性代数类的概述。这些对象呈现出两个接口：一个使对象看起来像一个具有全局索引操作的矩阵或向量，另一个使对象看起来像一个可以被单独处理的子块的集合。根据上下文，人们可能希望使用一个或另一个接口。



 通常，人们通过将构成物理量组的自由度（例如所有速度）归入线性系统的各个块来定义矩阵或向量的子结构。这在下面关于 @ref GlossBlock "块（有限元）"的词汇条中有更详细的定义。   </dd>







 <dt class="glossary">  @anchor  GlossBlock <b>Block (finite element)</b></dt>  <dd>  <i>Intent:</i> 块是 @ref GlossComponent  "组件 "的概括，它们将人们希望共同考虑的矢量值有限元的一个或多个组件组合起来。人们经常想这样做来定义与作用于向量值解的微分算子（部分）结构相对应的算子，例如  step-20  中的Schur补数求解器，或者  step-22  的块状求解器和预处理器。



 对于离散化的目的，块是更好的概念，因为并不总是能够解决一个解决方案的各个组成部分。特别是对于非  @ref  GlossPrimitive的 "原始 "元素来说，情况就是如此。以使用FE_RaviartThomas元素的混合拉普拉斯系统的解为例（见 step-20 ）。在那里，第一个<tt>dim</tt>分量是方向性速度。由于形状函数是这些的线性组合，这些<tt>dim</tt>分量只构成一个单一的块。另一方面，压力变量是标量，将形成第二个块，但在<tt>dim+1</tt>st分量中。



 每个块的最小尺寸是由底层有限元决定的（对于标量元素，一个块由一个分量组成，但以FE_RaviartThomas为例，一个块由<tt>dim</tt>分量组成）。然而，几个这样的最小块可以随意组合成用户定义的块，并根据应用情况进行组合。例如，对于<b>Q</b><sub>2</sub><sup><i>d</i></sup>-<b>Q</b><sub>1</sub>（Taylor-Hood）Stokes元素，有<i>d</i>+1个组件，原则上每个组件可以形成自己的块。但我们通常更感兴趣的是只有两个块，其中一个由所有的速度矢量分量组成（即这个块将有<i>d</i>分量），另一个只有单一的压力分量。



 <i>Implementation:</i> deal.II有许多不同的有限元类，它们都是从FiniteElement基类派生出来的（见 @ref feall  "有限元类模块"）。 除了一个例外，无论它们是标量还是矢量值，它们都定义了一个单一的块：有限元通过其 FiniteElement::n_components() 函数定义的所有矢量分量构成一个单一的块，即 FiniteElement::n_blocks() 返回一个。



 例外的是FESystem类，它采取多个较简单的元素，并将它们连接成较复杂的元素。因此，它可以有一个以上的块。一个FESystem有多少个块，就有多少个基础元素乘以它们的倍数（参见FESystem的构造函数来理解这一说法）。换句话说，它并不关心每个基元有多少个块，因此，你可以通过创建对象@code
    FESystem<dim> (FESystem<dim> (FE_Q<dim>(2), dim), 1,
                   FE_Q<dim>(1), 1);
 @endcode来产生一个只有两个块的斯托克斯元素。 另一方面，我们可以用@code
    FESystem<dim> (FE_Q<dim>(2), dim,
                   FE_Q<dim>(1), 1);
 @endcode来产生一个具有dim+1块的类似对象。
不过，除了块的数量之外，这两个对象在所有实际用途上都是一样的。



 <i>Global degrees of freedom:</i>
虽然我们在上面以矢量值解函数的矢量分量（或者，等同于以矢量值有限元空间）来定义块，但有限元的每个形状函数都是一个或另一个块的组成部分。因此，我们可以将定义在DoFHandler上的所有自由度划分为各个块。由于默认情况下DoFHandler类以一种或多或少的随机方式列举自由度，你首先要调用
DoFRenumbering::component_wise 函数以确保所有对应于单个块的自由度被连续列举。



 如果你这样做，你自然也会将矩阵和向量划分为块（见 @ref GlossBlockLA "块（线性代数）"）。 在大多数情况下，当你将矩阵或向量细分为块时，你是通过为有限元（即在大多数实际情况下的FESystem对象）定义的每个块创建一个块来实现的。然而，不必如此： DoFRenumbering::component_wise 函数允许将几个向量分量或有限元块归入同一个逻辑块（例如，见 @ref step_22 "  step-22 "或 step-31 教程程序，而不是 step-20 ）。因此，利用这一特性，我们可以实现同样的结果，即把矩阵细分为 $2\times 2$ 块，把向量细分为2块，对于上述创建斯托克斯元素的第二种方式，使用一个额外的参数，就像我们使用第一种创建有两个块的斯托克斯元素的方式一样。



 关于这个主题的更多信息可以在FESystem的文档中找到， @ref vector_valued 模块和其中参考的教程程序。



 <i>Selecting blocks:</i> 许多函数允许你将其操作限制在某些矢量分量或块上。例如，插值边界值的函数就是这种情况：人们可能只想插值有限元场的速度块的边界值，而不想插值压力块。这样做的方法是给这类函数传递一个BlockMask参数，见 @ref GlossBlockMask "本词汇表的block mask条目"。   </dd>







 <dt class=" glossary">  @anchor  GlossBlockMask <b>Block mask</b></dt>。



   <dd>  就像人们可以认为元素是由物理矢量分量（见 @ref GlossComponent ）或逻辑块（见 @ref GlossBlock ）组成一样，经常需要为不打算在有限元空间的<i>all</i>块上运行的操作选择一组此类块。使用BlockMask类来选择要操作的块。



 块掩码的工作方式与构件掩码基本相同，包括BlockMask类与ComponentMask类有类似的语义。更多信息请参见  @ref GlossComponentMask  "关于组件掩码的词汇表条目"。



   @note  虽然组件和块为具有多个向量分量的有限元提供了两种不同但同样有效的观点，但事实上，在整个库中，你可以传递ComponentMask参数而不是BlockMask参数的地方要多得多。幸运的是，一个可以转换为另一个，使用的语法 <code>fe.component_mask(block_mask)</code> where <code>block_mask</code> 是BlockMask类型的一个变量。换句话说，如果你有一个块掩码，但需要调用一个只接受组件掩码的函数，可以用这个语法来获得必要的组件掩码。



 <b>Creation of block masks:</b> 块掩码通常是通过要求有限元从某些选定的矢量分量中生成块掩码来创建的，使用这样的代码，我们创建了一个只表示斯托克斯元的速度分量的掩码（见 @ref vector_valued  ）。 @code
   FESystem<dim> stokes_fe (FESystem<dim>(FE_Q<dim>(2), dim), 1,    // Q2
element for the velocities FE_Q<dim>(1),                     1);     // Q1
element for the pressure FEValuesExtractors::Scalar pressure(dim); BlockMask
pressure_mask = stokes_fe.block_mask (pressure);
 @endcode 其结果是一个块状掩码，在1d以及2d和3d中，其值为  <code>[false,
true]</code>  。同样地，使用@code FEValuesExtractors::Vector velocities(0);
   BlockMask velocity_mask = stokes_fe.block_mask (velocities);
 @endcode在任何维度上都会产生一个掩码 <code>[true, false]</code> 。



 然而，请注意，如果我们以下列方式定义有限元。 @code
   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,    // Q2 element for the
velocities FE_Q<dim>(1), 1);     // Q1 element for the pressure
 @endcode，那么代码@code
   FEValuesExtractors::Scalar pressure(dim);
   BlockMask pressure_mask = stokes_fe.block_mask (pressure);
 @endcode将产生一个块状掩码，在2D中具有元素 <code>[false, false, true]</code> ，因为该元素具有 <code>dim+1</code> 成分和同样多的块。参见 @ref GlossBlock "本词汇表的块条目 "中关于块具体代表什么的讨论。   </dd>







 <dt class=" glossary">  @anchor  GlossBoundaryForm <b>Boundary form</b></dt>。



   <dd>
对于二维空间中的二维三角，边界形式是定义在面上的一个向量。它是单元格表面上坐标向量的图像的向量乘积。它是一个对面的法向矢量，指向外面，具有面元素的长度。



 一个更普遍的定义是（至少到这个矢量的长度为止），它正是考虑部分积分时必须的矢量，即形式为
$\int_\Omega \text{div} \vec \phi = -\int_{\partial\Omega} \vec n \cdot \vec
\phi$ 的等式。使用这个定义也解释了在嵌入空间 <code>spacedim</code> 的维数
<code>dim</code>
的域（和相应的三角形）的情况下，这个向量应该是什么：在这种情况下，边界形式仍然是一个定义在三角形面上的向量；它与边界的所有切线方向正交，并且在域的切线平面内。请注意，这与情况
<code>dim==spacedim</code> 是兼容的，因为那里的切平面是整个空间 ${\mathbb
R}^\text{dim}$  。



 在任何一种情况下，矢量的长度都等于参考面到当前单元面的变换行列式。   </dd>







 <dt class="glossary">  @anchor  GlossBoundaryIndicator <b>Boundary
indicator</b></dt>。



   <dd>  在Triangulation对象中，边界的每一部分都可以与一个唯一的数字（类型为
types::boundary_id)
，用于确定什么样的边界条件将被应用到边界的特定部分。边界是由单元格的面组成的，在三维中，是这些面的边缘。



 默认情况下，一个网格的所有边界指标都是零，除非你从一个网格文件中读取，并特别将其设置为不同的内容，或者你使用命名空间GridGenerator中的一个网格生成函数，该函数有 @ref GlossColorization "着色 "选项。一个典型的将部分边界指示器设置为其他东西的代码是这样的，这里将位于 $x=-1$ 的所有面的边界指示器设置为42：@code
   for (auto &face : triangulation.active_face_iterators())
     if (face->at_boundary())
       if (face->center()[0] == -1)
         face->set_boundary_id (42);
 @endcode 这调用函数 TriaAccessor::set_boundary_id.
在3D中，可能也适合在每个选定的面中调用 TriaAccessor::set_all_boundary_ids
。要查询某个特定面或边的边界指标，请使用 TriaAccessor::boundary_id. 。



 DoFTools和VectorTools命名空间中的许多函数都需要参数来指定边界的哪一部分，而且它们特别提到了boundary_ids。例如
DoFTools::make_periodicity_constraints,   DoFTools::extract_boundary_dofs,
DoFTools::make_zero_boundary_constraints  和
VectorTools::interpolate_boundary_values,
VectorTools::compute_no_normal_flux_constraints.  。



   @note  边界指标在网格细化时从母面和边继承到子面。关于边界指示器的更多信息，也在Triangulation类的文档中有所介绍。



   @note  对于 parallel::distributed::Triangulation, 类型的平行三角形，仅在开始时设置一次边界指标是不够的。参见  parallel::distributed::Triangulation  类的文档中关于这个主题的长篇讨论。   </dd>







 <dt class="glossary">  @anchor  GlossCoarseMesh <b>Coarse mesh</b></dt>  <dd>
deal.II中的 "粗略网格
"是一个仅由未被细化的单元组成的三角形对象，即一个没有单元是另一个单元的子单元的网格。这通常是deal.II中首次构建三角形的方式，例如，使用命名空间GridGenerator中的函数、GridIn类中的函数，或者直接使用函数
Triangulation::create_triangulation().
当然，我们可以在这样的网格上进行计算，但大多数时候（例如，参见几乎所有的教程程序），我们首先要全局地细化粗略的网格（使用
Triangulation::refine_global()),
]或自适应地细化（在这种情况下，首先计算一个细化准则，然后使用命名空间GridRefinement中的一个函数，最后调用
Triangulation::execute_coarsening_and_refinement()). ）。



   在某些情况下，我们也使用 "三角形的粗略网格 "这一短语，指的是三角形开始时的单元集，即所有目前 @ref GlossActive 三角形的 "活动单元 "都是通过网格细化得到的。当然，一些粗略的网格单元也可能是活跃的，如果它们从未被细化的话）。



   三角形对象以<i>levels</i>的方式存储单元：特别是，粗网格的所有单元都在零层。他们的子单元（如果我们在粗网格上执行
`Triangulation::refine_global(1)`
）将在第一层，等等。三角形的粗网格（在上一段的意义上）正好由三角形的零级单元组成。(它们是否处于活动状态(即没有子代)或已被细化，对这个定义并不重要)。



   deal.II中的大多数三角形类都存储了整个粗网格和至少一些细化的单元。( dealii::Triangulation 和 parallel::shared::Triangulation 类实际上都存储了整个网格的<i>all</i>个单元，而其他一些类如 parallel::distributed::Triangulation 在并行计算的每个进程中只存储了 @ref GlossActive 个 "活动单元 "的<i>some</i>个。) 在这些情况下，人们可以查询所有粗略网格单元的三角结构。其他三角剖分类（例如， parallel::fullydistributed::Triangulation) 只存储部分粗网格。这种情况也请参见 @ref GlossCoarseCellId "粗网格单元ID的概念"。   </dd>







 <dt class="glossary">  @anchor  GlossCoarseCellId <b>Coarse cell ID</b></dt>  <dd>  deal.II中的大多数三角形类，特别是 dealii::Triangulation,   parallel::shared::Triangulation,  和  parallel::distributed::Triangulation,  在并行计算的每个过程中存储三角形的整个 @ref GlossCoarseMesh  "粗略网格"。另一方面，其他类别则不是这样，特别是 parallel::fullydistributed::Triangulation, ，它是为粗略网格太大而不能存储在每个进程中，需要进行分区的情况而设计的。



   在这些情况下，在算法中经常需要唯一地引用一个粗略的网格单元。因为当前进程中的三角剖分对象实际上并没有存储整个粗网格，所以我们需要为每个粗网格单元设置一个全局唯一的标识符，这个标识符与本地存储的三角剖分零级中的索引无关。这个全局唯一的ID被称为
"粗略单元ID"。它可以通过函数调用@code
     triangulation.coarse_cell_index_to_coarse_cell_id (coarse_cell->index());
   @endcode进行访问，其中`triangulation`是指向零层单元的迭代器`coarse_cell`所属的三角结构。这里，`coarse_cell->index()`返回该单元在其细化层次中的索引（见
TriaAccessor::index()).
这是一个介于零和并行计算中当前进程上存储的粗网格单元数量之间的数字；它唯一地标识了该并行进程上的一个单元，但不同的并行进程可能对位于不同坐标的不同单元使用该索引。



   对于那些在每个进程上存储所有粗网格单元的类，
Triangulation::coarse_cell_index_to_coarse_cell_id()
只是返回可能的参数值的排列组合。在最简单的情况下，例如对于一个顺序的或并行的共享三角形，该函数实际上将简单地返回参数的值。对于其他情况，如
parallel::distributed::Triangulation,
，粗略单元ID的排序与粗略单元索引的排序不一样。最后，对于诸如
parallel::fullydistributed::Triangulation,
这样的类，该函数返回全局唯一的ID，它来自一个更大的可能指数集，而不是实际存储在当前进程上的粗放单元的指数。
</dd>







 <dt class="glossary">  @anchor  GlossColorization <b>Colorization</b></dt> <dd>
<em>  Colorization  </em>  是将三角图的某些部分用不同的标签标记的过程。颜色 <em>
一词的使用来自制图学，即通过给地图上的国家分配不同的颜色，使它们在视觉上相互区别。使用同一术语
<em> 着色 </em> 在数学中很常见，尽管我们给不同的区域分配整数而不是色调。



   <ol>   <li>  GridGenerator命名空间中的大多数函数都采取一个可选的参数  <code>colorize</code>  。这个参数控制边界的不同部分是否会被分配不同的  @ref GlossBoundaryIndicator  "边界指标"。   一些函数也会分配不同的  @ref GlossMaterialId  "材料指标"。 </li>   <li>  函数 GraphColoring::make_graph_coloring() 计算一个三角形的分解（更确切地说，是一个迭代器的范围）。没有两个相邻的单元被赋予相同的颜色。 </li>   </ol>   </dd> 。







 <dt class=" glossary">  @anchor  GlossComponent <b>Component</b></dt>



   <dd>  当考虑方程组时，其中的解不仅仅是一个单一的标量函数，我们说我们有一个<i>vector system</i>与一个<i>vector-valued solution</i>。例如，在 step-8 中考虑的弹性方程的矢量解是 $u=(u_x,u_y,u_z)^T$ ，由三个坐标方向上的位移组成。然后，该解决方案有三个元素。同样， step-22 中考虑的三维斯托克斯方程有四个元素。  $u=(v_x,v_y,v_z,p)^T$  . 我们在交易二中称矢量值解的元素为<i>components</i>。为了得到良好的解决，对于解有 $n$ 个元素，需要有 $n$ 个偏微分方程来描述它们。这个概念在 @ref vector_valued 模块中有详细的讨论。



 在有限元程序中，人们经常想解决这个矢量值解决方案的单个元素（组件），或组件的集合。例如，我们在 step-8 中做了大量的工作，在 @ref vector_valued "处理矢量值问题 "模块中也提供了大量的文档。如果你只考虑偏微分方程（不考虑其离散化），那么<i>components</i>的概念是自然的。



 另一方面，当谈论有限元和自由度时，<i>components</i>并不总是正确的概念，因为组件并不总是可以单独处理。特别是 @ref GlossPrimitive 的 "非原始有限元 "就是这种情况。同样，人们也不一定总是<i>want</i>处理单个分量，而是处理分量的集合&mdash；例如，所有速度分量在一起，并与斯托克斯系统中的压力分开，而不进一步将速度分割成其单个分量。在这两种情况下，正确的概念是 @ref GlossBlock "块 "的概念。 由于每个分量，如果可以单独处理，也是一个块，所以用块的方式思考是最常见的更好的策略。



 对于一个给定的有限元，可以使用 FiniteElementData::n_components()
函数查询组件的数量，可以使用 FiniteElement::get_nonzero_components().
找出对于一个给定的有限元形状函数，哪些向量组件是非零的。
形状函数的各个组件的值和梯度（如果该元素是原始的）可以使用参考单元上的
FiniteElement::shape_value_component() 和 FiniteElement::shape_grad_component()
函数查询。 FEValues::shape_value_component() 和 FEValues::shape_grad_component()
函数在实数单元上做同样的事情。也请参见FiniteElement和FEValues类的文档。



 <i>Selecting components:</i> 许多函数允许你将其操作限制在某些向量组件或块上。例如，插值边界值的函数就是这种情况：人们可能只想插值有限元场的速度分量的边界值，而不想插值压力分量。这样做的方法是给这类函数传递一个ComponentMask参数，见 @ref GlossComponentMask "本词汇表的组件掩码条目"。   </dd>







 <dt class=" glossary">  @anchor  GlossComponentMask <b>Component mask</b></dt>



   <dd>  当使用矢量值元素（见 @ref vector_valued ）来解决方程组时，人们经常希望将一些操作限制在只有某些解决变量。例如，在求解斯托克斯方程时，人们可能希望只插值速度分量的边界值而不插值压力。在deal.II中，这通常是通过传递函数a<i>component mask</i>完成的。分量掩码总是被指定为ComponentMask对象，我们可以把它看作一个数组，其条目数与有限元的分量一样多（例如，在Stokes情况下，有 <code>dim+1</code> 个分量），每个条目要么为真，要么为假。 在这个例子中，我们只想插值斯托克斯系统的速度分量的边界值，那么这个分量掩码将是 <code>[true, true, false]</code> in 2d and <code>[true, true, true, false]</code> 的3d，表示不应设置压力变量的边界值（解决方案中 <code>dim+1</code> 的最后一个矢量分量。



 有许多函数采取这样的分量掩码，例如 DoFTools::make_zero_boundary_values,
VectorTools::interpolate_boundary_values,   KellyErrorEstimator::estimate,
等。在某些情况下，有多个具有这些名称的函数，但只有其中一些具有分量掩码参数。



 <b>Semantics of component
masks:</b>许多函数，接受一个已经默认构建的分量掩码对象，表示<i>all
components</i>，即，就像向量有正确的长度，并且只填充了 <code>true</code>
的值。原因是默认初始化的对象可以使用代码片断 <code>ComponentMask()</code>
来构建到位，因此可以作为函数签名中的默认参数。



 换句话说，ComponentMask对象可以处于两种状态中的一种。它们可以被一个非零长度的布尔运算向量初始化；在这种情况下，它们代表一个特定长度的掩码，其中一些元素可能是真，另一些可能是假。或者，ComponentMask可能已经被默认初始化了（使用默认构造函数），在这种情况下，它代表了一个长度不确定的数组（即适合情况的长度），其中<i>every
entry</i>为真。



 <b>Creation of component masks:</b> 分量掩码通常是通过要求有限元从某些选定的分量中生成一个分量掩码来创建的，使用这样的代码，我们创建一个掩码，只表示斯托克斯元的速度分量（见 @ref vector_valued  ）。 @code
   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,    // Q2 element for the
velocities FE_Q<dim>(1), 1);     // Q1 element for the pressure
   FEValuesExtractors::Scalar pressure(dim);
   ComponentMask pressure_mask = stokes_fe.component_mask (pressure);
 @endcode 其结果是一个分量掩码，在2d中，它的值是 <code>[false, false,
true]</code>  。同样地，使用@code FEValuesExtractors::Vector velocities(0);
   ComponentMask velocity_mask = stokes_fe.component_mask (velocities);
 @endcode会导致2d中的掩码 <code>[true, true, false]</code>
。当然，在3D中，结果将是 <code>[true, true, true, false]</code>  。



   @note  正如我们可以认为组成元素是由 @ref GlossComponent "组件 "或 @ref GlossBlock "块 "组成的一样，有组件掩码（由ComponentMask类表示）和 @ref GlossBlockMask "块掩码"（由BlockMask类表示）。 FiniteElement类有在这两种对象之间转换的函数。



   @note  并非所有的构件掩码实际上都有意义。例如，如果你有一个2d中的FE_RaviartThomas对象，那么有一个 <code>[true, false]</code> 形式的分量掩码是没有任何意义的，因为你试图选择一个有限元的单个矢量分量，其中每个形状函数都有 $x$ 和 $y$ 速度。从本质上讲，虽然你当然可以创建这样的分量掩码，但你对它无能为力。   </dd>











 <dt class=" glossary">  @anchor  GlossCompress <b>Compressing distributed
                                              vectors and matrices</b></dt>。



   <dd>
对于%并行计算，deal.II使用PETScWrappers和TrilinosWrappers命名空间中定义的向量和矩阵类。当使用MPI在%parallel中运行程序时，这些类只在当前处理器上存储一定数量的行或元素，而向量或矩阵的其余部分则存储在属于我们MPI宇宙的其他处理器上。当你组装线性系统时，这就出现了一定的问题：我们向矩阵和右手边的向量添加元素，这些元素可能在本地存储，也可能不在。有时，我们也可能只想<i>set</i>一个元素，而不是向它添加。



 PETSc和Trilinos都允许添加或设置没有本地存储的元素。在这种情况下，它们将我们想要存储或添加的值写入缓存，我们需要调用其中一个函数
TrilinosWrappers::VectorBase::compress(),
TrilinosWrappers::SparseMatrix::compress(),
PETScWrappers::VectorBase::compress()  或  PETScWrappers::MatrixBase::compress()
，然后将缓存中的值运送到拥有应该被添加或写入的元素的MPI进程。由于MPI模型只允许从发送方发起通信（也就是说，它不是一个远程过程调用），这些函数是集体的，也就是说，它们需要被所有处理器调用。



 然而，有一个障碍：PETSc和Trilinos都需要知道这些 <code>compress()</code> 函数所调用的操作是适用于添加元素还是设置元素。 在某些情况下，并不是所有的处理器都会添加元素，例如，当使用非常 @ref GlossCoarseMesh 的 "粗略（初始）网格 "时，一个处理器不拥有任何单元。 出于这个原因，compress()需要一个VectorOperation类型的参数，它可以是::%add，或者::%insert。 从7.3版本开始，这个参数对于向量和矩阵是必需的。



 简而言之，你需要在以下情况下调用compress()（而且只在这些情况下，虽然在其他情况下调用compress()只是花费一些性能）。



 1. 在你对矩阵和向量的汇编循环结束时。如果你直接写条目或使用
AffineConstraints::distribute_local_to_global.  使用 VectorOperation::add.
，就需要这样做。



 2.
当你完成了对矩阵/向量中单个元素的设置，然后再进行其他操作（向元素添加，其他操作如缩放、求解、读取等）。使用
VectorOperation::insert. 。



 3. 像2.一样，但用于向单个元素加值。使用 VectorOperation::add. 。



 所有其他的操作，如缩放或添加向量，赋值，调用deal.II（VectorTools,
AffineConstraints, ...）或求解器不需要调用compress()。   </dd>



   @note  压缩是一个只适用于向量的操作，其元素在一个并行的MPI宇宙中被一个且唯一的处理器拥有。它不适用于 @ref GlossGhostedVector "有幽灵元素的向量"。







 <dt class=" glossary">  @anchor  GlossConcept <b>Concepts in deal.II</b></dt>



   <dd>  在deal.II中，有几个地方我们要求模板中的类型与某个接口相匹配或以某种方式表现：这种约束在C++中被称为 <em>  概念 </em> 。更多信息请参见 @ref Concepts 中的讨论和deal.II中的概念列表。   </dd>







 <dt class=" glossary">  @anchor  GlossDimension <b>Dimensions `dim` and
`spacedim`</b></dt>。



   <dd>  deal.II中的许多类和函数都有两个模板参数， @p dim 和 @p spacedim.  一个例子是基本的Triangulation类。 @code
   template <int dim, int spacedim=dim>
   class Triangulation {...};
 @endcode
在所有这些上下文中，当你看到`dim`和`spacedim`被引用时，这些参数有如下含义。



   <ul>   <li>   @p dim  表示网格的维度。例如，一个由线段组成的网格是一维的，因此对应于`dim==1`。由四边形组成的网格为`dim==2`，六面体的网格为`dim==3`。 </li>



     <li>   @p spacedim  表示这种网格所在空间的维度。一般来说，一维网格生活在一维空间中，同样，二维和三维网格也是如此，它们将二维和三维领域细分。因此， @p  spacedim模板参数的默认值等于 @p dim. ，但情况并不一定如此。例如，我们可能想解决地球表面的沉积物迁移方程。在这种情况下，域是地球的二维表面（`dim==2`），它生活在三维坐标系中（`spacedim==3`）。 </li>   </ul>  。



 更一般地说，deal.II可以用来解决嵌入高维空间的<a
href="https://en.wikipedia.org/wiki/Manifold">manifolds</a>上的偏微分方程。换句话说，这两个模板参数需要满足`dim
<= spacedim'，尽管在许多应用中，我们只需满足`dim == spacedim'。



 按照几何学的惯例，我们说 "二维
"被定义为`spacedim-dim`。换句话说，一个由四边形组成的三角形，其坐标是三维的（我们将使用`Triangulation<2,3>`对象）具有
"二维"。



 这两个参数不一样的使用例子在  step-34  ,  step-38  ,  step-54  中显示。   </dd>







 <dt class="glossary">  @anchor  GlossDoF <b>Degree of freedom</b></dt>。



   <dd>  术语 "自由度"（通常缩写为
"DoF"）在有限元界通常用来表示两个略有不同但相关的事情。首先是我们想把有限元解表示为形状函数的线性组合，形式为
$u_h(\mathbf{x}) = \sum_{j=0}^{N-1} U_j \varphi_j(\mathbf{x})$  。这里， $U_j$
是一个膨胀系数的向量。因为我们还不知道它们的值（我们将计算它们作为线性或非线性系统的解），它们被称为
"未知数 "或
"自由度"。该术语的第二个含义可以解释如下。对有限元问题的数学描述通常是说，我们正在寻找一个满足某些方程组的有限维函数
$u_h \in V_h$ （例如， $a(u_h,\varphi_h)=(f,\varphi_h)$ 的所有测试函数
$\varphi_h\in V_h$ ）。换句话说，我们在这里说的是，解决方案需要位于某个空间
$V_h$
中。然而，为了在计算机上实际解决这个问题，我们需要选择这个空间的一个基；这就是我们在上面用系数
$U_j$ 展开 $u_h(\mathbf x)$ 时使用的形状函数 $\varphi_j(\mathbf{x})$
的集合。当然，空间 $V_h$
的基数有很多，但我们将特别选择由传统上在网格单元上局部定义的有限元函数描述的基数。在这种情况下描述
"自由度 "需要我们简单地 <i>enumerate</i> 空间的基函数  $V_h$  。对于 $Q_1$
元素，这意味着简单地以某种方式列举网格的顶点，但对于更高的元素，还必须列举与网格的边、面或单元内部有关的形状函数。提供这种列举
$V_h$ 的基础函数的类被称为DoFHandler。 列举自由度的过程在deal.II中被称为
"分配DoF"。   </dd>



 <dt class=" glossary">  @anchor  GlossDirectionFlag <b>Direction
flags</b></dt>。



   <dd>  <i>direction
flag</i>用于嵌入高维空间的三角形中，以表示单元的方向，并使流形具有导向性。它使用
CellAccessor::direction_flag()
进行访问，并在创建三角化时由三角化类进行设置。你可以使用
Triangulation::flip_all_direction_flags() 函数来改变一个三角形的所有方向标志。



 这个标志对于像这样的情况是必须的：假设我们有一个嵌入二维空间的一维网格。



     @image html direction_flag.png "One dimensional mesh in two dimensions"



 在一维空间的一维网格中，我们总是可以确保一个单元的左边顶点的位置比右边顶点的位置的值要小。然而，如果我们将网格嵌入到高维空间中，我们就不能再这样做了。例如，上面的网格中的单元格可以用下面的顶点集来描述。<code>(0,1),
(1,2), (3,2), (4,3), (4,5)
</code>。(作为附带说明，注意这里我们有一些顶点--例如顶点2--是不止一个单元的右端点。)
如果我们将每个单元的法线定义为与连接线的第一个顶点和第二个顶点的矢量垂直的单位矢量，那么我们最终会得到以下图片。



     @image html direction_flag_normals.png "Normal vectors"



 换句话说，这个一维流形是没有方向的。我们原则上可以在创建这样的网格时恢复顶点的顺序（尽管有很好的理由不这样做，例如，这个网格可能是由提取二维网格的表面网格产生的，而我们希望保留每个线段的顶点顺序，因为它们目前与二维单元的面的顶点顺序一致）。在deal.II中选择的另一种策略是简单地与每个单元关联，法线应该是该单元的左法线还是右法线。在上面的例子中，五个单元格的标志将是<code>true,
true, false, false,
true</code>。根据每个单元上的标志值，将右法线与正负1相乘，可以得到一组为流形定位的法向量。



 类似的问题发生在三维空间的二维网格上。我们注意到，如果二维流形不是可定向的，就不可能找到一致的方向标志；这样的流形目前不被deal.II支持。
</dd>







 <dt class=" glossary">  @anchor  GlossDistorted <b>Distorted cells</b></dt>



   <dd>  <i>distorted
cell</i>是指从参考单元到实数单元的映射有一个雅各布系数，其行列式在单元的某处为非正值的单元。通常情况下，我们只在单元格的顶点检查这个行列式的符号。函数
GeometryInfo::alternating_form_at_vertices 可以计算这些顶点的行列式。



 举例来说，如果所有的行列式都是大致相等的值，并且在 $h^\text{dim}$
的顺序上，那么这个单元格就是好的形状。例如，一个正方形单元或面的行列式等于
$h^\text{dim}$
，而一个强剪切的平行四边形的行列式则小得多。同样地，一个边长很不相等的单元格会有差异很大的行列式。反之，一个被夹住的单元，其中两个或更多顶点的位置被折叠成一个点，这个位置的行列式为零。最后，一个倒置或扭曲的单元，其中两个顶点的位置是失序的，将有负的行列式。



 下面两张图片显示了2D和3D的一个完好单元、一个捏合单元和一个扭曲的单元。



   @image html distorted_2d.png "A well-formed, a pinched, and a twisted cell in 2d."



   @image html distorted_3d.png "A well-formed, a pinched, and a twisted cell in 3d."



 扭曲的单元格可以以两种不同的方式出现。原始的 @ref GlossCoarseMesh "粗略网格 "可能已经包含了这样的单元，或者它们可能是由于移动或扭曲了一个相对较大的网格而产生的结果。



 如果在创建三角网格时给出适当的标志，由GridGenerator和GridIn中的各种函数调用的函数
Triangulation::create_triangulation, （但也可以由用户代码调用，见 step-14 和
step-49 末尾的例子 ]，将通过抛出一个类型为 Triangulation::DistortedCellList.
的异常来提示创建带有变形单元的粗略网格。
如果你不打算在这些单元上装配任何东西，创建带有变形单元（特别是塌陷/针状单元）的网格是合法的。例如，考虑这样一种情况：人们想模拟一种有液体填充的裂缝的弹性材料的行为，如一个储油罐。如果压力过大，裂缝就会闭合--而离散裂缝体积的单元就会塌陷成零体积。只要你不在这些单元上进行积分来模拟流体的行为（如果裂缝的体积为零，就不存在任何流体），这种网格是完全合法的。因此，
Triangulation::create_triangulation
不是简单地中止程序，而是抛出一个包含变形单元列表的异常；这个异常可以被捕获，如果你认为你可以忽略这个条件，你可以通过对捕获的异常不做任何反应。



 在某些情况下，函数 GridTools::fix_up_distorted_child_cells
可以通过移动具有未变形父单元的变形子单元的顶点，来修复精化网格上的变形单元。



 请注意，Triangulation类在默认情况下不会测试是否存在扭曲的单元，因为确定一个单元是否扭曲并不是一个便宜的操作。如果你想让Triangulation对象测试单元格的变形，你需要在创建对象时通过传递适当的标志来指定这一点。
</dd>







 <dt class="glossary">  @anchor  distributed_paper <b>Distributed computing
paper</b></dt>。



   <dd>  "分布式计算论文 "是W. Bangerth, C. Burstedde, T. Heister和M.
Kronbichler的一篇论文，题为
"大规模并行通用有限元代码的算法和数据结构"，描述了deal.II中%并行分布式计算的实现，即不仅像
step-17
中的线性系统被分割到不同机器的计算，还包括三角计算和DoFHandler对象。实质上，它是
parallel::distributed  命名空间和  step-40  中所用技术的指南。



 该论文的完整参考资料如下。 @code{.bib}
@Article{BBHK11,
  author =       {Wolfgang Bangerth and Carsten Burstedde and Timo Heister
                  and Martin Kronbichler},
  title =        {Algorithms and data structures for massively parallel generic
  adaptive finite element codes},
  journal =      {ACM Trans. Math. Softw.},
  year =         2011,
  volume =       38,
  pages =        {14/1--28}}
 @endcode



 对于大规模的%并行计算，deal.II建立在<a href="http://www.p4est.org/"
target="_top">p4est</a>库的基础上。如果你使用这个功能，也请引用他们网站上列出的p4est论文。
</dd>







 <dt class="glossary">  @anchor  GlossFaceOrientation <b>Face
orientation</b></dt>  <dd>  在一个三角形中，通过应用右手边的规则（x,y ->
normal），可以从面的方向推导出法向量。
我们注意到，在2D的标准方向中，面0和面2的法线指向单元格，面1和面3的法线指向外部。在3D中，面0、2和4的法线指向单元格内，而面1、3和5的法线指向外面。这些信息同样可以从
GeometryInfo<dim>::unit_normal_orientation. 中查询到。



 然而，事实证明，相当数量的三维网格不能满足这个约定。这是由于一个单元的面的约定已经暗示了相邻单元的东西，因为它们共享一个共同的面，对第一个单元的固定也固定了两个单元的相对面的法向量。很容易构建单元格循环的案例，对于这些案例，我们无法为所有的面找到与这个约定一致的方向。



 由于这个原因，上述惯例只是我们所说的 <em> 标准方向 </em>  ......II实际上允许3d中的面具有标准方向，或者其相反的方向，在这种情况下，构成单元格的线会有还原的顺序，法向量会有相反的方向。你可以通过调用<tt>cell->face_orientation(face_no)</tt>来询问一个单元是否有标准方向：如果结果是 @p true, ，那么这个面有标准方向，否则它的法向量就会指向另一个方向。在应用程序中，你需要这些信息的地方其实并不多，但库中有几个地方用到了这个。请注意，在2D中，结果总是 @p true. 。然而，虽然2D中的每个面总是在标准方向上，但有时你可以指定一些东西来假设不是这样的；一个例子是函数 DoFTools::make_periodicity_constraints(). 。



 还有两个描述面的方向的标志：face_flip和face_rotation。在GeometryInfo类中存在一些关于这些的文档。
DoFTools::make_periodicity_constraints 函数中给出了它们在用户代码中的使用实例。
</dd>







 <dt class="glossary">  @anchor  GlossGeneralizedSupport <b>Generalized support points</b></dt>  <dd>  "Generalized support points"，顾名思义，是 @ref GlossSupport  "support points "的泛化。后者是用来描述一个有限元在个别点（"支持点"）上的简单<i>interpolates</i>值。如果我们称这些点为 $\hat{\mathbf{x}}_i$ （其中帽子表示这些点定义在参考单元上 $\hat{K}$ ），那么人们通常以这样的方式定义形状函数 $\varphi_j(\mathbf{x})$ ，即<i>nodal functionals</i> $\Psi_i[\cdot]$ 只是在支持点评估函数，即 $\Psi_i[\varphi]=\varphi(\hat{\mathbf{x}}_i)$ ，并且选择基础以便 $\Psi_i[\varphi_j]=\delta_{ij}$  其中 $\delta_{ij}$  是Kronecker三角函数。这就导致了常见的 @ref GlossLagrange "拉格朗日元素"。



 (在矢量值的情况下，除了支持点 $\hat{\mathbf{x}}_i$
外，唯一需要提供的其他信息是<i>vector component</i>  $c(i)$
第1个节点函数对应的，因此 $\Psi_i[\varphi]=\varphi(\hat{\mathbf{x}}_i)_{c(i)}$
。)



 另一方面，还有其他种类的元素不是这样定义的。例如，对于最低阶的Raviart-Thomas元素（见FE_RaviartThomas类），节点函数评估的不是一个具有 @p dim 分量的矢量值有限元函数的各个分量，而是这个矢量的<i>normal component</i>。   $\Psi_i[\varphi]
  =
  \varphi(\hat{\mathbf{x}}_i) \cdot \mathbf{n}_i
 $  ，其中 $\mathbf{n}_i$ 是 $\hat{\mathbf{x}}_i$
所在的单元格面的法向量。换句话说，当在 $\hat{\mathbf{x}}_i$ 处评估时，节点函数是
$\varphi$ 组件的<i>linear
combination</i>。类似的事情也发生在BDM、ABF和Nedelec元素上（见FE_BDM、FE_ABF、FE_Nedelec类）。



 在这些情况下，元素没有 <i>support points</i>
，因为它不是纯粹的插值；然而，在定义形状函数时，仍然涉及某种插值，因为节点函数仍然需要在特殊点上进行点评估
$\hat{\mathbf{x}}_i$  。 在这些情况下，我们称这些点为<i>generalized support
points</i>。



 最后，还有一些元素仍然不适合这个方案。例如，一些层次化的基函数（例如，见FE_Q_Hierarchical元素）的定义是这样的：节点函数是有限元函数的<i>moments</i>，2d的
$\Psi_i[\varphi]
  =
  \int_{\hat{K}} \varphi(\hat{\mathbf{x}})
  {\hat{x}_1}^{p_1(i)}
  {\hat{x}_2}^{p_2(i)}
 $ ，同样，3d的 $p_d(i)$ 是形状函数 $i$
描述的矩的顺序。其他一些元素使用边或面的矩。在所有这些情况下，节点函数根本不是通过插值定义的，那么这些元素既没有支持点，也没有广义支持点。
</dd>







 <dt class="glossary">  @anchor  geometry_paper <b>geometry paper</b></dt>  <dd>
"geometry paper "是L. Heltai, W. Bangerth, M. Kronbichler, and A.
Mola的一篇论文，题目是 "Using exact geometry information in finite element
computations"，它描述了
deal.II如何描述域的几何结构。特别是，它讨论了Manifold类所基于的算法基础，以及它需要为网格细化、法向量的计算和其他许多几何学进入有限元计算的地方提供什么样的信息。



 这篇论文目前可在arXiv网站https://arxiv.org/abs/1910.09824。
这篇论文的完整参考资料如下。 @code{.bib}
@misc{heltai2019using,
    title={Using exact geometry information in finite element computations},
    author={Luca Heltai and Wolfgang Bangerth and Martin Kronbichler and Andrea
Mola}, year={2019}, eprint={1910.09824}, archivePrefix={arXiv},
    primaryClass={math.NA}
}
 @endcode  </dd>  。







 <dt class="glossary">  @anchor  GlossGhostCell <b>Ghost cells</b></dt>  <dd>  如果一个网格使用 parallel::distributed::Triangulation 类分布在多个MPI进程中，每个处理器只存储自己拥有的单元，其他处理器拥有的一层相邻单元，所有 @ref GlossCoarseMesh  "粗级单元"，以及所有为维持相邻单元必须最多相差一个细化级别这一不变性所需的单元。存储在每个进程上的不属于该进程但与该进程拥有的单元相邻的单元被称为 "幽灵单元"，对于这些单元，谓词 <code>cell-@>is_ghost()</code> 返回真。幽灵单元被保证存在于全局分布的网格中，也就是说，这些单元实际上被另一个进程所拥有，并且在那里没有被进一步细化。更多信息请参见 @ref distributed_paper  "分布式计算论文"。



 幽灵单元层包括所有与任何本地拥有的单元相邻的面、边或顶点的单元，而这些单元本身并不是本地拥有的。换句话说，幽灵单元完全包围了本地拥有的单元的子域（当然，域的边界除外）。



 幽灵单元的概念对于在每个处理器上存储整个网格的三角计算没有意义，即三角计算和
parallel::shared::Triangulation 类。   </dd>







 <dt class="glossary">  @anchor  GlossGhostedVector <b>Ghosted vectors</b></dt>
<dd>
在并行计算中，向量一般有两种：没有和有鬼元素。没有鬼魂元素的向量在处理器之间唯一地划分了向量元素：每个向量条目正好有一个处理器拥有它，而且这个处理器是唯一存储这个条目的值的。换句话说，如果零号处理器存储了一个向量的0...49号元素，一号处理器存储了50...99号元素，那么一号处理器访问这个向量的42号元素就不走运了：它没有被存储在这里，也无法评估其值。
这将导致一个断言。



 另一方面，在很多情况下，我们需要知道不属于本地的向量元素，例如在本地拥有的单元上评估解决方案（见 @ref GlossLocallyOwnedCell ），其中一个自由度位于我们不属于本地的单元的界面上（在这种情况下，它必须是 @ref GlossGhostCell  ]"幽灵单元"），而相邻的单元可能是其所有者--换句话说，该自由度不是 @ref GlossLocallyOwnedDof "本地拥有的"，而只是 @ref GlossLocallyActiveDof "本地活动的DoFs"。这种自由度的值通常存储在拥有该自由度的机器上，因此，在当前机器上无法访问。



 因为人们经常需要这些值，所以有第二种矢量，通常称为
"幽灵矢量"。幽灵向量在每个处理器上存储一些元素，而该处理器不是所有者。
对于这样的向量，你可以读取你当前所在的处理器所存储的那些元素，但你不能写入这些元素，因为要做到这一点，需要将新的值传播给所有其他拥有这个值副本的处理器（这些处理器的列表可能是当前处理器不知道的，也没有办法有效地找到）。因为你不能写进重影向量，所以初始化这样一个向量的唯一方法是通过从一个非重影向量的赋值。这意味着我们必须从其他处理器中导入那些我们本地想要存储的元素。



 幽灵向量的实际存储方式在并行向量的各种实现中是不同的。对于PETSc（以及相应的
PETScWrappers::MPI::Vector
类），重影向量存储的元素与非重影向量相同，另外还有一些由其他处理器拥有的额外元素。换句话说，对于每个元素，在所有的处理器中都有一个明确的所有者，那些当前处理器存储但不拥有的元素（即
"幽灵元素"）只是其他地方的主值的镜像--因此，被称为 "幽灵"。
parallel::distributed::Vector 类也是这种情况。



 另一方面，在Trilinos中（以及因此在 TrilinosWrappers::MPI::Vector),
中，鬼魂向量仅仅是元素分布重叠的平行向量的一个视图。幽灵化的
"Trilinos向量本身不知道哪些条目是幽灵化的，哪些是局部拥有的。事实上，一个重影向量甚至可能不会存储所有非重影向量在当前处理器上会存储的元素。因此，对于Trilinos向量来说，不存在我们在非鬼魂情况下（或在PETSc情况下）所拥有的向量元素的
"所有者 "的概念，"鬼魂元素
"这个名字可能会产生误导，因为在这个观点中，我们在本地可用的每个元素可能也会被存储在其他地方，但即使是这样，本地元素也不是主位置的镜像值，因为每个元素没有所有者。



   @note   @ref distributed  文档模块简要介绍了不同种类的向量通常用于何处。   </dd>







 <dt class="glossary">  @anchor  hp_paper <b>%hp-paper</b></dt>  <dd>  "hp-paper
"是W. Bangerth和O. Kayser-Herold的一篇论文，题目是
"hp有限元软件的数据结构和要求"，它描述了在实现deal.II的hp-framework时所使用的许多算法和数据结构。特别是，它总结了许多使用连续元素的%hp-有限元必须考虑的棘手问题。



 这篇论文的完整参考资料如下。 @code{.bib}
@Article{BK07,
  author =       {Wolfgang Bangerth and Oliver Kayser-Herold},
  title =        {Data Structures and Requirements for hp Finite Element
                  Software},
  journal =      {ACM Trans. Math. Softw.},
  year =         2009,
  volume =       36,
  number =       1,
  pages =        {4/1--4/31}
}
 @endcode可从<a
href="http://www.math.colostate.edu/~bangerth/publications.html">http://www.math.colostate.edu/~bangerth/publications.html</a>获得，也可参见<a
href="https://www.dealii.org/publications.html#details">deal.II
publications</a>了解详情。



 该论文中显示的数值例子是用稍加修改的  step-27
版本生成的。与该教程程序的主要区别是，该程序中的各种操作都是为论文计时的，以比较不同的选项，并表明
$hp$ 方法真的不是那么昂贵。   </dd>







 <dt class="glossary">  @anchor  GlossInterpolation <b>Interpolation with finite elements</b></dt>  <dd>  用有限元插值的目的是计算代表有限元函数的系数向量，从而使原始函数和有限元函数的 @ref GlossNodes  "节点值 "相吻合。因此，插值过程包括对给定函数<i>f</i>的所有 @ref GlossNodes "节点函数"<i>N<sub>i</sub></i>进行评估，并将结果作为条目<i>i</i>存储在系数向量中。   </dd>







 <dt class="glossary">  @anchor  GlossLagrange <b>Lagrange elements</b></dt>  <dd>  基于拉格朗日插值的有限元在  @ref GlossSupport  "支持点"。   </dd>







 <dt class="glossary">  @anchor  GlossLocallyOwnedCell <b>Locally owned cell</b></dt>  <dd>  当使用分布式网格时，这个概念标识了所有单元的子集，见 @ref distributed 模块。在这样的网格中，每个单元正好被一个处理器所拥有。本地拥有的是那些由当前处理器拥有的。



 并行计算中的每个处理器都有一个覆盖整个域的三角形，它由本地拥有的单元、 @ref GlossGhostCell "幽灵单元 "和 @ref GlossArtificialCell "人工单元 "组成。   </dd>







 <dt class="glossary">  @anchor  GlossLocallyOwnedDof <b>Locally owned degrees of freedom</b></dt>  <dd>  当使用分布式网格时，这个概念标识了所有自由度的一个子集，见 @ref distributed 模块。 本地拥有的自由度生活在本地拥有的单元上。由于自由度只属于一个处理器，不同处理器所拥有的单元之间的接口上的自由度可能属于一个或另一个处理器，所以并非本地拥有的单元上的所有自由度也是本地拥有的自由度。



 本地拥有的自由度是 @ref GlossLocallyActiveDof "本地活动自由度 "的一个子集。   </dd>







 <dt class="glossary">  @anchor  GlossLocallyActiveDof <b>Locally active degrees of freedom</b></dt>  <dd>  当使用分布式网格时，这个概念标识了所有自由度的一个子集，见 @ref distributed  模块。 本地活动的自由度是那些生活在本地拥有的单元上的自由度。因此，在不同处理器所拥有的单元之间的界面上的自由度属于一个以上处理器的本地活动自由度集合。



 本地活动自由度是 @ref GlossLocallyOwnedDof "本地拥有的自由度 "的超集和 @ref GlossLocallyRelevantDof "本地相关自由度 "的子集。   </dd>







 <dt class="glossary">  @anchor  GlossLocallyRelevantDof <b>Locally relevant degrees of freedom</b></dt>  <dd>  这个概念在使用分布式网格时确定了所有自由度的子集，见 @ref distributed  模块。 本地相关的自由度是那些生活在本地拥有的或幽灵单元上的自由度。 因此，它们可能被不同的处理器所拥有。



 本地相关自由度是 @ref GlossLocallyActiveDof "本地活动自由度 "的一个超集。   </dd>







 <dt class=" glossary">  @anchor  GlossManifoldIndicator <b>%Manifold
indicator</b></dt>。



   <dd>  构成三角网格的每个对象（单元、面、边等），都与一个唯一的编号（类型为
types::manifold_id) ，用于识别哪个流形对象在网格细化时负责生成新的点。



 默认情况下，一个网格的所有流形指标都被设置为 numbers::flat_manifold_id.
。一个典型的设置对象上的流形指标的代码会是这样的，这里为所有中心的 $x$
分量小于0的单元设置流形指标为42。



 @code
 for (auto &cell : triangulation.active_cell_iterators())
   if (cell->center()[0] < 0)
     cell->set_manifold_id (42);
 @endcode



 这里我们调用函数 TriaAccessor::set_manifold_id(). ，也可以调用
TriaAccessor::set_all_manifold_ids
，以递归地设置每个面（如果是三维面，还有边）的流形指标。要查询某个特定对象边缘的流形指标，请使用
TriaAccessor::manifold_id(). 。



 上面的代码只是设置了Triangulation的特定部分的流形指标，但它本身并没有改变Triangulation类在网格细化中对待这个对象的方式。为此，你需要调用
Triangulation::set_manifold()
来将流形对象与特定的流形指标联系起来。这允许Triangulation对象使用不同的方法在单元格、面或边上寻找新的点来进行细化；默认情况是对所有面和边使用FlatManifold对象。



   @note  网格细化时，流形指标会从父辈继承到子辈。关于流形指示器的更多信息，也在Triangulation类的文档部分以及 @ref manifold  "流形文档模块 "中介绍。歧管指标在  step-53  和  step-54  中使用。   </dd>



   @see   @ref manifold  "关于歧管的模块"







 <dt class="glossary">  @anchor  GlossMaterialId <b>Material id</b></dt>  <dd>
三角形的每个单元都有一个名为 "材料ID
"的属性与之相关。它通常用于具有异质系数的问题，以确定一个单元在域的哪一部分，因此，系数应该在这个特定的单元上具有哪个值。在实践中，一个单元的材料ID通常用于识别哪些单元属于域的特定部分，例如，当你有不同的材料（钢铁、混凝土、木材），但都属于同一个域。在组装双线性表格的过程中，我们通常会查询与某个单元相关的材料ID，并使用它来确定（例如，通过表格查询，或一连串的if-else语句）该单元的正确材料系数是多少。



 这个材料ID可以在构建三角形时设置（通过CellData数据结构），也可以在之后通过使用单元格迭代器设置。关于这个功能的典型使用，请看
step-28
的教程程序。GridGenerator命名空间的函数通常将所有单元的材料ID设置为0。当通过GridIn类读取三角图时，不同的输入文件格式有不同的约定，但通常是明确指定材料ID，如果没有，则GridIn简单地将其设置为零。因为一个单元的材料是与域的一个特定区域相关的，材料ID在网格细化时由子单元从父单元继承。



 材料ID的设置和查询使用 CellAccessor::material_id, CellAccessor::set_material_id
和  CellAccessor::recursively_set_material_id  函数。   </dd>







 <dt class="glossary">  @anchor  GlossMPICommunicator <b>MPI Communicator</b></dt>  <dd>  在消息传递接口（MPI）的语言中，一个通信器可以被认为是一个邮件系统，允许向邮件系统的其他成员发送消息。在每个通信器中，每个 @ref GlossMPIProcess "进程 "都有一个 @ref GlossMPIRank "等级"（相当于门牌号），允许识别消息的发送方和接收方。不可能通过通信器向不属于该通信器/邮件服务的接收者发送消息。



 当通过诸如@code
  mpirun -np 32 ./step-17
 @endcode（或在你的集群上使用的批处理提交系统中使用的等价物）的命令行调用启动一个并行程序时，MPI系统启动32份
step-17 可执行文件。 其中每个都可以访问 <code>MPI_COMM_WORLD</code>
通信器，然后由所有32个处理器组成，每个都有自己的等级。这个MPI宇宙中的一个进程子集后来可以同意创建其他通信器，只允许在一个进程子集之间进行通信。
</dd>







 <dt class="glossary">  @anchor  GlossMPIProcess <b>MPI Process</b></dt>  <dd>
在分布式内存机器上运行并行作业时，人们几乎总是使用MPI。在那里，像@code mpirun
-np 32 ./step-17
 @endcode这样的命令行调用（或者在你的集群上使用的批处理提交系统中使用的等价物）启动32份
step-17
可执行文件。其中一些实际上可能在同一台机器上运行，但一般来说，它们将在不同的机器上运行，不能直接访问对方的内存空间。



 在消息传递接口（MPI）的语言中，在（可能是不同的）机器上运行的同一可执行文件的每个副本被称为<i>processes</i>。所有并行运行的进程的集合被称为 "MPI宇宙"，并由  @ref GlossMPICommunicator  "MPI通信器"  <code>MPI_COMM_WORLD</code>  来识别。



 每个进程只能立即访问其自身内存空间中的对象。一个进程不能从其他进程的内存中读取或写入。因此，进程可以通信的唯一方式是互相发送消息。也就是说（正如在
step-17
的介绍中所解释的那样），人们通常会调用更高级别的MPI函数，而作为通信器一部分的所有进程都参与其中。一个例子是计算一组整数的总和，每个进程提供总和的一个项。
</dd>







 <dt class="glossary">  @anchor  GlossMPIRank <b>MPI Rank</b></dt>  <dd>  在消息传递接口（MPI）的语言中， @ref GlossMPIProcess  "MPI进程 "的<i>rank</i>是该进程在当前作为一个并行作业运行的所有进程的集合 <code>MPI_COMM_WORLD</code> 中所携带的数字。更正确地说，它是 @ref GlossMPICommunicator "MPI通信器 "中的数字，该通信器将所有进程的一个子集与一个并行作业组合在一起（其中 <code>MPI_COMM_WORLD</code> 简单地表示<i>complete</i>进程集）。



 在每个通信器中，每个进程都有一个独特的等级，与所有其他进程的等级不同，可以在MPI通信调用中识别一个接收方或发送方。在一个处理器上运行的每个进程可以通过调用
Utilities::MPI::this_mpi_process(). 来查询自己在通信器中的等级。
参与通信器的进程总数（即通信器的<i>size</i>）可以通过调用
Utilities::MPI::n_mpi_processes().  </dd> 获得。







 <dt class="glossary">  @anchor  mg_paper <b>%Multigrid paper</b></dt>  <dd>
"multigrid paper "是B. Janssen和G. Kanschat的一篇论文，题目是 "Adaptive
Multilevel Methods with Local Smoothing for H1- and Hcurl-Conforming High Order
Finite Element Methods"，描述了在实现
deal.II的多网格框架时所用的许多算法和数据结构。它是实现 step-16
中用于多网格方法的类的基础。



 本文的完整参考资料如下。 @code{.bib}
@article{janssen2011adaptive,
  title=    {Adaptive Multilevel Methods with Local Smoothing for H^1- and
H^{curl}-Conforming High Order Finite Element Methods}, author=   {Janssen,
B{\"a}rbel and Kanschat, Guido}, journal=  {SIAM Journal on Scientific
Computing}, volume=   {33}, number=   {4}, pages=    {2095--2114}, year= {2011},
  publisher={SIAM}}
 @endcode 论文见<a
href="http://dx.doi.org/10.1137/090778523">DOI:10.1137/090778523</a>，更多细节见<a
href="https://www.dealii.org/publications.html#details">deal.II
publications</a>。   </dd>







 <dt class=" glossary">  @anchor  GlossNodes <b>Node values or node
functionals</b></dt>。



   <dd>  习惯上将有限元定义为三要素 $(K,P,\Psi)$  其中 -  $K$  是单元，在deal.II中这总是一个线段、四边形或六面体； -  $P$  ]是一个有限维空间，例如，从 @ref GlossReferenceCell "参考单元 "映射到 $K$ 的多项式空间； -  $\Psi$ 是一组 "节点函数"，即函数 $\Psi_i : P \rightarrow {\mathbb R}$  。  $P$ 的维度必须等于节点函数的数量。 有了这个定义，我们可以定义局部函数空间的基础，即一组 "形状函数" $\varphi_j\in P$ ，要求 $\Psi_i(\varphi_j) = \delta_{ij}$ ，其中 $\delta$ 是克朗克三角。



 这种对有限元的定义有几个优点，涉及分析和实施。对于分析来说，它意味着与某些空间
(FiniteElementData::Conformity),
的一致性，例如连续性，是由节点函数决定的。在deal.II中，它有助于大大简化像FE_RaviartThomas这样的复杂元素的实现。



 节点函数的例子是 @ref GlossSupport "支持点 "的值和关于Legendre多项式的矩。例子。



   <table><tr>
   <th>Element</th>
   <th>%Function space</th>
   <th>Node values</th></tr>
   <tr><th>FE_Q, FE_DGQ</th>
     <td>[1.x.109]</td>
     <td>values in support points</td></tr>
   <tr><th>FE_DGP</th>
     <td>[1.x.110]</td>
     <td>moments with respect to Legendre polynomials</td></tr>
   <tr><th>FE_RaviartThomas (2d)</th>
     <td>[1.x.111]</td>
     <td>moments on edges and in the interior</td></tr>
   <tr><th>FE_RaviartThomasNodal</th>
     <td>[1.x.112]</td>
     <td>Gauss points on edges(faces) and anisotropic Gauss points in the
interior</td></tr>
 </table>



 如上所述，有限元的构造允许编写描述有限元的代码，只需提供一个多项式空间（不需要给它任何特定的基础--只要方便就完全足够）和节点函数。例如，在
FiniteElement::convert_generalized_support_point_values_to_dof_values()
函数中就用到了这一点。   </dd>







 <dt class="glossary">  @anchor  GlossParallelScaling <b>Parallel
scaling</b></dt>  <dd>  当我们说一个并行程序可以 "扩展
"时，我们的意思是，如果我们让它解决的问题变大，程序不会变得过于缓慢（或占用过多的内存），如果我们保持问题大小不变但增加处理它的处理器（或核）数量，运行时间和内存消耗将按比例下降。



 更具体地说，想想一个问题，其大小由一个数字 $N$
给出（可以是单元格的数量，未知数的数量，或其他一些指示性的数量，如解决它所需的CPU周期的数量），对于这个问题，你有
$P$ 个处理器可用于解决。在一个理想的世界里，这个程序需要的运行时间是 ${\cal
O}(N/P)$
，这意味着我们可以通过提供更多的处理器将运行时间减少到任何想要的值。同样，为了使程序具有可扩展性，其总体内存消耗需要为
${\cal O}(N)$ ，在每个参与的进程上需要为 ${\cal O}(N/P)$
，这再次意味着我们可以通过提供足够多的处理器，将任何问题纳入计算机附加在每个处理器上的固定内存量。



 对于可扩展性的实际评估，我们经常区分 "强 "和 "弱
"可扩展性。这些评估渐进式的声明，如 ${\cal O}(N/P)$ 极限运行时间 $N\rightarrow
\infty$ 和/或 $P\rightarrow \infty$  。具体来说，当我们说一个程序是 "强可扩展性
"时，我们的意思是，如果我们有一个固定大小的问题 $N$
，那么我们可以通过向该问题投掷更多的处理器来减少运行时间和内存消耗（在每个处理器上）与
$P$
成反比。特别是，强可扩展性意味着，如果我们提供两倍的处理器，那么每个进程的运行时间和内存消耗都将减少2倍。换句话说，通过提供越来越多的处理器，我们可以越来越快地解决<i>same
problem</i>的问题。



 相反，"弱可扩展性 "意味着，如果我们将问题大小 $N$
增加一个固定的系数，并将可用于解决问题的处理器 $P$
的数量增加相同的系数，那么总体运行时间（以及每个处理器的内存消耗）保持不变。换句话说，我们可以通过提供越来越多的处理器，在相同的壁时钟时间内解决<i>larger
and larger problems</i>。



 在这个理论意义上，没有程序是真正可扩展的。相反，一旦 $N$ 或 $P$
的增长超过一定的限度，所有的程序就不再具有可扩展性。 因此，我们经常说
"程序可以扩展到4000个核心"，或者 "程序可以扩展到100,000,000个未知数
"这样的话。程序不能无限制扩展的原因有很多；只要看看（相对简单的） step-17
的教程程序，就可以说明这些原因了。 -
顺序性的部分。许多程序都有不能或不能并行化的代码部分，也就是说，一个处理器必须做一定的、固定的工作量，不会因为周围有总共
$P$ 个处理器而减少。在 step-17
中，生成图形输出时就是这种情况：一个处理器为整个问题创建图形输出，也就是说，它需要做
${\cal O}(N)$ 工作。这意味着这个函数的运行时间为 ${\cal O}(N)$ ，而不考虑 $P$
，因此整个程序将无法达到 ${\cal O}(N/P)$ 的运行时间，而是有一个可以描述为
$c_1N/P + c_2N$
的运行时间，其中第一项来自可扩展的操作，如组装线性系统，而后者来自在进程0上生成图形输出。如果
$c_2$
足够小，那么程序可能看起来对小数量的处理器具有强扩展性，但最终强扩展性将停止。此外，程序也不能弱扩展，因为增加问题的大小
$N$ ，同时以同样的速度增加处理器的数量 $P$
，并不能保持这一个函数的运行时间不变。 - 重复的数据结构。在 step-17
中，每个处理器存储整个网格。也就是说，每个处理器都要存储一个大小为 ${\cal O}(N)$
的数据结构，不管 $P$
如何。最终，如果我们使问题的大小足够大，即使我们增加处理器的数量，这也会溢出每个处理器的内存空间。因此，很明显，这种复制的数据结构可以防止程序弱速扩展。
但它也阻止了程序的强扩展，因为为了创建一个大小为 ${\cal O}(N)$
的对象，至少要写到 ${\cal O}(N)$ 的内存位置，要花费 ${\cal O}(N)$
的CPU时间。因此，如果我们提供越来越多的处理器，整个算法的一个组成部分就不会表现得
${\cal O}(N/P)$ 。 - 通信。仅举一个例子，如果你想计算一个向量的 $l_2$
常数，所有MPI进程都存储了一些条目，那么每个进程都需要计算其自身条目的平方和（这将需要
${\cal O}(N/P)$
时间，因此可以完美扩展），但随后每个进程都需要将其部分和发送到一个进程，将它们全部相加并取平方根。在最好的情况下，发送一个包含单个数字的信息需要恒定的时间，而不考虑进程的总体数量。因此，同样地，每一个做通信的程序都不能强势扩展，因为程序中有些部分的CPU时间要求并不随着你为固定规模分配的处理器数量而减少
$P$
。在现实中，情况实际上更糟糕：参与一个通信步骤的进程越多，一般来说需要的时间就越长，例如，因为要把每个人的贡献加起来的那个进程必须把所有的东西加起来，需要
${\cal O}(P)$
时间。换句话说，CPU的时间<i>increases</i>与进程的数量有关，因此不仅阻止了程序的强扩展，而且也阻止了弱扩展。实际上，MPI库并不通过将每个消息发送到一个进程，然后将所有的东西加起来来实现
$l_2$ 规范；相反，它们在树上做成对的减少，而不是像 ${\cal O}(P)$
那样增长运行时间，而是像 ${\cal O}(\log_2 P)$
那样，以发送更多消息为代价。尽管如此，根本的一点是，当你增加更多的处理器时，运行时间将以
$P$ 的方式增长，而不管操作的实际实现方式如何，因此它不能扩展。)



 这些以及其他阻碍程序完美扩展的原因可以在<a
href="https://en.wikipedia.org/wiki/Amdahl%27s_law"> <i>Amdahl's
law</i><i>Amdahl's law</i></a>中总结出来，即如果一个程序的整体工作 $\alpha$
的一部分可以被并行化，即可以在 ${\cal O}(\alpha W/P)$
时间内运行，而该程序工作的一部分 $1-\alpha$
不能被并行化（即。它包括只有一个进程可以做的工作，例如在 step-17
中生成图形输出；或者每个进程都必须以复制的方式执行，例如将带有本地贡献的消息发送到一个专门的进程进行积累），那么程序的总体运行时间将是@f{align*}
   T = {\cal O}\left(\alpha \frac WP + (1-\alpha)W \right).
 @f} ] 因此，你得到的 "加速"，即你的程序在 $P$
处理器上的运行速度与在单个进程上运行程序相比的系数（假设这是可能的）将是@f{align*}
   S = \frac{W}{\alpha \frac WP + (1-\alpha)W}
     = \frac{P}{\alpha + (1-\alpha)P}.
 @f} 如果 $\alpha<1$ ，它对所有实际存在的程序都是如此，那么 $S\rightarrow
\frac{1}{1-\alpha}$ 就是 $P\rightarrow \infty$
，这意味着有一个点，在这个问题上投入更多处理器不会再有任何显著回报了。



 在实践中，重要的是<i>up to which problem size</i>或<i>up to which number of
processes</i>或<i>down to which size of local problems
${\cal}(N/P)$</i>一个程序的规模。对于deal.II，经验表明，在大多数具有合理快速网络的集群上，人们可以解决多达几十亿个未知数的问题，最多有几千个处理器，而每个进程的未知数则在40,000到100,000之间。最后一个数字是最相关的：如果你有一个问题，例如
$10^8$
未知数，那么在1000-2500个处理器上解决它是有意义的，因为每个进程处理的自由度数量保持在40000以上。因此，每个进程都有足够的工作要做，所以
${\cal O}(1)$
的通信时间并不占优势。但是用1万个或10万个处理器来解决这样的问题是没有意义的，因为这些处理器的每个局部问题都变得非常小，以至于它们的大部分时间都是在等待通信，而不是在做自己部分的工作。
</dd>



 <dt class="glossary">  @anchor  GlossPeriodicConstraints <b>Periodic boundary
 conditions</b></dt>  <dd>
周期性边界条件经常在只有部分物理相关域被建模时使用。人们假设解决方案只是在被认为是周期性的边界上周期性地继续。
在deal.II中，通过 DoFTools::make_periodicity_constraints() 和
GridTools::collect_periodic_faces(). 支持这一点。一旦使用
parallel::distributed::Triangulation ，还必须调用
parallel::distributed::Triangulation::add_periodicity()
以确保所有进程知道周期性边界两边的三角化的相关部分。一个典型的分布式三角剖分的过程是。
-# 创建一个网格 -# 使用 GridTools::collect_periodic_faces()
（三角形）收集周期性面 -# 使用
parallel::distributed::Triangulation::add_periodicity() 将周期性信息添加到网格中
-# 使用 GridTools::collect_periodic_faces() （DoFHandler）收集周期性面 -# 使用
DoFTools::make_periodicity_constraints() 添加周期性约束。



 这方面的一个例子可以在  step-45  中找到。   </dd>







 <dt class="glossary">  @anchor  GlossPrimitive <b>Primitive finite
 elements</b></dt>  <dd>  如果从形状函数数到矢量 @ref  GlossComponent "分量 "之间有唯一的关系，那么有限元（由其形状函数描述）就是原始的。这意味着，如果一个元素是原始的，那么矢量值元素的每个形状函数正好有一个非零分量。这尤其包括所有标量元素以及通过FESystem类从其他基元（例如标量）元素组装的矢量值元素，如 step-8 、 step-29 、 step-22 和其他一些元素所示。另一方面， step-20 和 step-21 中使用的FE_RaviartThomas类或FE_Nedelec类提供了非原始有限元，因为在那里，每个矢量值形状函数可能有几个非零分量。   </dd>







 <dt class="glossary">  @anchor  GlossReferenceCell <b>Reference cell</b></dt>
<dd>
超立方体[0,1]<sup>dim</sup>，所有参数化的有限元形状函数都在其上定义。参考单元的许多属性由GeometryInfo类描述。
</dd>







 <dt class="glossary">  @anchor  GlossSerialization <b>Serialization</b></dt>。



   <dd>  术语 "序列化
"指的是将一个对象的状态写入流中，之后再检索的过程。一个典型的用例是将程序的状态保存到磁盘上，以便以后可能的复活，通常是在长期运行的计算的检查点/重启策略的背景下，或者在不是很可靠的计算机上（例如，在非常大的集群上，个别节点偶尔会出现故障，然后导致整个MPI作业的崩溃）。在这两种情况下，人们希望偶尔保存程序的状态，以便在失败时，可以在该点重新启动，而不是从头开始运行。



 deal.II通过实现<a
 href="http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/index.html"
 target="_top">BOOST
serialization</a>库的必要接口实现了序列化设施。关于如何保存和恢复对象的例子见那里。
</dd>







 <dt class="glossary">  @anchor  GlossShape <b>Shape functions</b></dt>  <dd>
将有限元基函数限制在单个网格单元。   </dd>







 <dt class="glossary">  @anchor  GlossSubdomainId <b>Subdomain id</b></dt>  <dd>  三角形的每个单元都有一个称为 "subdomain id "的属性，可以通过像  <code>cell-@>subdomain_id()</code>  这样的调用进行查询，例如，可以通过  <code>cell-@>set_subdomain_id(13)</code>  进行设置。(这些调用分别解析为 CellAccessor::subdomain_id() 和 CellAccessor::set_subdomain_id(), 。)虽然原则上这个属性可以被应用程序认为有用的任何方式使用（它只是一个与每个单元相关的整数，可以表示任何你想要的），至少对于以%并行方式运行的程序，它通常表示 "拥有 "这个单元的处理器的 @ref GlossMPIRank "MPI等级"。



 对于基于MPI并行化的程序，但每个处理器都存储整个三角形（例如， step-17 和
step-18 ，但不是 step-40 ），子域ID通过划分网格分配给单元，然后每个MPI进程只对它
"拥有
"的单元工作，即。属于处理器拥有的子域（传统上，这是子域ID的情况，其数值与MPI通信器中MPI进程的等级一致）。分区通常使用
GridTools::partition() 函数完成，但也可以使用任何其他方法来完成。(另外，
parallel::shared::Triangulation 类可以使用类似的方法自动划分网格)。



 另一方面，对于使用MPI并行化的程序，但网格是使用 parallel::distributed::Triangulation 类在多个处理器上保持分布的，单元的子域ID与拥有该单元的处理器相联系。换句话说，查询一个单元的子域id可以告诉你这个单元是否被当前的处理器所拥有（即如果<code>cell-  @>subdomain_id()  ==  triangulation.parallel::distributed::Triangulation::locally_owned_subdomain()</code>)  或者被其他处理器所拥有。在%并行分布的情况下，子域ID只分配给当前处理器拥有的单元以及紧邻的 @ref GlossGhostCell "幽灵单元"。 更远的单元被保留在每个处理器上，以确保每个MPI进程都能访问完整的 @ref GlossCoarseMesh "粗网格"，并确保相邻的单元最多相差一个细化级别这一不变性。这些单元被称为 "人工"（见 @ref GlossArtificialCell "这里"），具有特殊的子域id值 types::artificial_subdomain_id. 。



 除了常规的子域id外，还有第二套密切相关的标志，与每个单元相关。"水平子域id"。
这些标志不仅存在于活动单元，事实上，也存在于网格层次结构中的每一个单元。它们的含义完全类似于常规的子域id，但它们是由
CellAccessor::level_subdomain_id() 和 CellAccessor::set_level_subdomain_id()
函数读写的。   </dd>







 <dt class="glossary">  @anchor  GlossSupport <b>Support points</b></dt>  <dd>
支持点根据定义是那些  $p_i$  ，使得对于形状函数  $v_j$  持有  $v_j(p_i) =
\delta_{ij}$  。因此，有限元插值可以由支持点中的值唯一地定义。



 拉格朗日元素填充由 FiniteElement::get_unit_support_points(), 访问的矢量，从而使函数 FiniteElement::has_support_points() 返回<tt>真</tt>。当然，这些支持点是在 @ref GlossReferenceCell 的 "参考单元 "上。 然后，可以使用FEValues（与Mapping结合使用）来访问实际网格单元上的支持点。



   @note   @ref GlossSupport  "支持点 "的概念只限于基于拉格朗日插值的有限元族。对于一个更普遍的概念，见 @ref GlossGeneralizedSupport "广义支持点"。   </dd>







 <dt class="glossary">  @anchor  GlossTargetComponent <b>Target
component</b></dt>  <dd>
当向量和矩阵按分量分组为块时，往往希望将几个原始分量集合为一个分量。例如，这可能是将斯托克斯系统的速度分组为一个单一的块。
</dd>







 <dt class="glossary">  @anchor  GlossUnitCell <b>Unit cell</b></dt>  <dd>  见 @ref GlossReferenceCell  "参考单元"。   </dd>







 <dt class="glossary">  @anchor  GlossUnitSupport <b>Unit support points</b></dt>  <dd>  这些是参考单元上的 @ref GlossSupport "支持点"，定义在FiniteElement中。例如，1d中通常的Q1元素在<tt>x=0</tt>和<tt>x=1</tt>有支持点（类似地，在更高维度的单位正方形或立方体的顶点上）。另一方面，高阶拉格朗日元素的单位支持点也在单位线、正方形或立方体的内部。   </dd>







 <dt class="glossary">  @anchor  GlossUserFlags <b>User flags</b></dt>  <dd>
三角形为每条线、四边形等提供一个用户标志位。
这个字段可以像所有其他数据一样使用迭代器进行访问，使用的语法是@code
      cell->set_user_flag();                // set the user flag of a cell
      if (cell->user_flag_set() == false)   // if cell hasn't been flagged yet
        {
           cell->face(0)->set_user_flag();  // flag its first face
        }
   @endcode。
通常，如果一个算法走过所有的单元，并且需要另一个单元，例如邻居，是否已经被处理过的信息，那么这个用户标志就可以使用。同样，它也可以用来标记边界上的面、四边形或线，对它们已经进行了一些操作。后者通常是有用的，因为像@code
      // in 3d
      for (cell=dof_handler.begin_active();
           cell!=dof_handler.end(); ++cell)
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          if (cell->line(l)->at_boundary())
            {
               do something with this line
            }
   @endcode这样的循环会不止一次地遇到一些边界线。因此，人们会在循环的主体中设置该线的用户标志，并且只有在用户标志先前没有被设置的情况下才会进入主体。有许多额外的函数可以通过迭代器接口访问；更多信息请参见TriaAccessor类。请注意，没有用户标志可以与顶点相关联；然而，由于顶点是连续编号的，这可以很容易地在用户代码中使用一个bools矢量来模拟。



   有两个函数， Triangulation::save_user_flags 和 Triangulation::load_user_flags
，可以从一个流或一个bools向量中写入和读取这些标志。与
Triangulation::save_refine_flags 和 Triangulation::load_refine_flags,
不同的是，这两个函数存储和读取所有使用过的线、四边形等的标志，也就是说，不仅仅是活动的标志。



   如果你想存储更多具体的用户标志，你可以使用函数
Triangulation::save_user_flags_line 和 Triangulation::load_user_flags_line
，对四边形等也是如此。



   至于细化和粗化标志，这些函数有两个版本，一个是从流中读/写，一个是从<tt>向量中读/写
@<bool@></tt>. ，后者用于临时存储标志，而第一个用于将它们存储在文件中。



   在使用前用 Triangulation::clear_user_flags()
函数清除用户标志是一个很好的做法，因为经常需要在多个函数中使用这些标志。如果在调用需要这些标志的函数时，这些标志可能还在使用中，那么这个函数应该按照上述方法保存和恢复这些标志。



     @note  如果需要在单元格、行或面中存储更多的信息，而不仅仅是一个布尔标志，那么请参见 @ref GlossUserData  "用户数据"。   </dd>







 <dt class="glossary">  @anchor  GlossUserData <b>User pointers and user indices</b></dt>  <dd>  就像 @ref GlossUserFlags "用户标志 "一样，Triangulation类为每条线、四边形和六边形提供了一个字段，在其中存储比单一布尔标志更具描述性的数据。这被称为 "用户数据"，可以在其中存储的数据是一个无符号整数或一个空指针。两者都是用来索引到一个更大的数组中，该数组包含了应用程序想要附加到网格实体的更详细的数据。



   用户数据的存储和检索方式如下。   @code
      for (cell=dof_handler.begin_active();
           cell!=dof_handler.end(); ++cell)
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          if (cell->line(l)->at_boundary())
            {
              cell->line(l)->set_user_index(42);
            }
   @endcode 同样地，有函数 TriaAccessor::set_user_pointer 用来设置指针，
TriaAccessor::user_index 和 TriaAccessor::user_pointer
用来检索索引和指针。要清除所有的用户索引或指针，请使用
Triangulation::clear_user_data().
与标志一样，有一些函数允许保存和恢复用户数据，可以是网格层次结构的所有实体，也可以是线、四边形或六边形的单独数据。还有一些额外的函数可以通过迭代器接口访问；更多信息请参见TriaAccessor类。



     @note  用户指针和用户索引被存储在同一个地方。为了避免不必要的转换，Triangulation会检查它们中的哪一个正在使用，并且不允许访问另一个，直到 Triangulation::clear_user_data() 被调用。



     @note  关于 @p void 指针的类型安全缺失的常规警告显然在这里得到了体现；类型的正确性等的责任完全在于指针的使用者。   </dd>







 <dt class=" glossary">  @anchor  workstream_paper <b>%WorkStream paper</b></dt>
<dd>  "WorkStream paper "是一篇由 B. Turcksin、M. Kronbichler 和 W. Bangerth
写的论文，讨论了 WorkStream
的设计和实现。WorkStream的核心是一种设计模式，即在有限元代码中反复使用的东西，因此，可以通用地实现。特别是，本文阐述了这种模式的动机，然后提出了实现它的不同方法。它还比较了不同实现方式的性能。



 这篇论文的完整参考资料如下。 @code{.bib}
@Article{TKB16,
  author =       {Bruno Turcksin and Martin Kronbichler and Wolfgang Bangerth},
  title =        {\textit{WorkStream} -- a design pattern for multicore-enabled
finite element computations}, journal =      {accepted for publication in the
ACM Trans. Math. Softw.}, year =         2016
}
 @endcode 可从<a
href="http://www.math.colostate.edu/~bangerth/publications.html">http://www.math.colostate.edu/~bangerth/publications.html</a>中获取，详情还可参见<a
href="https://www.dealii.org/publications.html#details">deal.II
publications</a>。   </dd>







 <dt class="glossary">  @anchor  GlossZOrder <b>Z order</b></dt>  <dd>  单元的
"Z顺序 "描述了一个单元被遍历的顺序。



  默认情况下，如果你在deal.II中的所有单元上写一个循环，单元将按照较粗的单元（即从 @ref GlossCoarseMesh "粗网格 "单元中以较少的细化步骤获得的单元）在较细的单元（即以较多的细化步骤获得的单元）之前的顺序进行遍历。在每个细化层中，单元的遍历顺序与它们的创建顺序有关；然而，从本质上讲，这个顺序最好被认为是 "未指定的"：你将以某种顺序精确访问一个给定细化层中的每个单元，但你不应该对这个顺序做出任何假设。



  因为单元格的创建顺序会影响到单元格的顺序，所以对于两个相同的网格，你遍历单元格的顺序可能会发生变化。例如，想想一个有两个单元的1d（粗）网格。如果你先精炼其中的第一个单元，然后再精炼另一个单元，那么你将以不同的顺序遍历精炼层1上的四个单元，而不是先精炼第二个粗单元，再精炼第一个粗单元。



  这种顺序对于几乎所有的应用都是完全实用的，因为在大多数情况下，以何种顺序遍历单元实际上并不重要。此外，它允许使用导致特别低的高速缓存缺失频率的数据结构，因此对高性能计算应用来说是有效的。



  另一方面，在某些情况下，人们希望以特定的、可重复的顺序遍历单元，这只取决于网格本身，而不是其创建历史或任何其他看似任意的设计决定。Z顺序
"是实现这一目标的方法之一。



  为了解释Z顺序的概念，考虑下面的网格序列（每个单元使用 "级别.索引 "符号进行编号，其中 "级别 "是指从 @ref GlossCoarseMesh "粗略网格 "单元到特定单元所需的细化次数，"索引 "是该单元在特定细化级别中的索引）。



    @image html simple-mesh-0.png "A coarse mesh"   @image html simple-mesh-1.png "The mesh after one refinement cycle"   @image html simple-mesh-2.png "The mesh after two refinement cycles"   @image html simple-mesh-3.png "The mesh after three refinement cycles" 。



  注意第2层的单元格是如何按照它们被创建的顺序排列的。这并不总是这样：如果中间有单元格被移除，那么新创建的单元格就会填上这样产生的洞）。



  那么deal.II遍历单元格的 "自然 "顺序就是0.0->1.0->1.1->1.2->1.3->2.0->2.1->2.2->2.3->2.4->2.5->2.6->2.7（如果你只想遍历 @ref GlossActive "活动单元格"，那么省略这个列表中所有有子单元。  这可以被认为是对数字 "level.index "的 "lexicographic "顺序，但由于每个级别内的索引没有被很好地定义，这不是一个特别有用的概念。另外，我们也可以把它看作是对应于这个网格的树的一个可能的广度优先遍历，它代表了单元格之间的父子关系。



    @image html simple-mesh-tree.png "The tree that corresponds to the mesh after three refinement cycles"



  另一方面，Z顺序对应于树的一个特定的深度优先的遍历。即：从一个单元格开始，如果它有孩子，那么就遍历这些单元格的孩子；只要一个孩子有孩子，这个规则就会递归应用。



  对于上述给定的网格，这产生了以下顺序。0.0 -> 1.0 -> 2.4 -> 2.5 -> 2.6 -> 2.7
-> 1.1 -> 1.2 -> 1.3 -> 1.4 -> 2.0 -> 2.1 -> 2.2
-> 2.3（同样，如果你只关心活动单元，那么从这个列表中删除0.0，1.0，和1.3）。因为一个单元格的子代顺序是明确定义的（相对于每一层内的单元格顺序），这种
"分层 "遍历是有意义的，尤其是独立于三角测量的历史。



  在实践中，使用一个递归函数可以很容易地实现。  @code
    template <int dim>
    void visit_cells_hierarchically (const typename
Triangulation<dim>::cell_iterator &cell)
    {
      if (cell->has_children())
        for (unsigned int c=0; c<cell->n_children(); ++c)
          visit_cells_hierarchically (cell->child(c));
      else
        {
          ... do whatever you wanted to do on each cell ...;
        }
    }
  @endcode 然后按以下方式调用这个函数。  @code
    // loop over all coarse mesh cells
    for (typename Triangulation<dim>::cell_iterator cell =
triangulation.begin(0); cell != triangulation.end(); ++cell)
      visit_cells_hierarchically (cell);
  @endcode



  最后，作为对 "Z
"顺序一词的解释：如果你按照这种分层方式出现的顺序画一条穿过所有单元的线，那么它在每个细化的单元上看起来就像一个左-右倒置的Z。事实上，这样定义的曲线可以被认为是一条空间填充曲线，有时也被称为
"莫顿排序"，见https://en.wikipedia.org/wiki/Z-order_curve 。   </dd>











   </dl>  */
