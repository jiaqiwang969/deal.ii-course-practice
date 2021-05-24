# deal.ii-course-practice
这是个人的课后练习，是对Theory and Practice of Finite Element Methods这门课程的巩固和强化。 
- 课程链接<https://www.youtube.com/playlist?list=PLcvf2raG3YsHis9dJfQ-anJgR1b7xGum3> 
- 课程导师 Luca Heltai  <luca.heltai@sissa.it>
- 代码文档: <http://deal-ii.com/deal.ii-course-practice/>



## 练习 01 - 了解代码的组成结构

目标：
- 建立一个现代化的编程环境
- 熟悉git、VSCode、CMake、GoogleTest和C++。

任务：
- 将源文件`pythagoras.cc`分成两个文件，一个包含`main`函数，另一个包含
`main`函数，另一个包含测试。
- 把它们移到`source`目录下，并命名为`main.cc`和`point_tests.cc`。
- 修改`CMakeLists.txt`文件，以确保程序仍能编译和运行
- 增加三个测试，以检查三维版的`Point<dim>`类是否按预期运行


## 练习 02 - Finite Elements的理论与实践 
目标：
- 认识和了解Triangulation, DoFHandler, FiniteElement的基础类

任务：
- 你将需要以下程序：
-- 安装了 "remote development"插件的Visual Studio Code；
-- Paraview (版本 > 5.4 会更好)；
-- 检查文件`.devcontainer`目录是否包含适合你的配置。当你打开版本库的根目录时，它应该会提示你在容器中重新打开。如果没有，请尝试使用 "F1 -> Reopen folder in container "命令进行配置；
-- 如果在容器中打开，你应该能够运行终端(`CTRL + ~`)并看到类似的东西。

``` 
要以管理员（用户 "root"）身份运行一个命令，使用 "sudo <command>"。
详情见 "man sudo_root"。
dealii@3f25621373bd:/workspaces/02_tria_dh_fem$ 
```
- 该文件夹包含一个`CMakeLists.txt`文件。VSCode应该检测到这一点，并在终端窗口的底部提供配置和构建应用程序的命令。


- 如果你有ubuntu/debian并且不想使用Docker的话，可以进行如下设置。
 -- 在这种情况下，你可以复制docker镜像中的命令，你可以在这里找到：
 <https://github.com/dealii/docker-files/blob/master/dependencies-focal/Dockerfile>
来安装所有的依赖项。
--在这之后，你可以使用下面的命令安装最新版本的的库，使用这里的命令： <https://github.com/dealii/dealii/blob/master/contrib/docker/Dockerfile>

### Lab-02: step-1

1.  阅读步骤1的文档，网址是 <https://www.dealii.org/current/doxygen/deal.II/step_1.html>

2. 在VSCode容器内编译并运行`step-1`，并查看输出。
3. 了解如果你现在一个接一个地调用`first_grid(triangulation)`和`second_grid(triangulation)`会怎样。修复你得到的问题，不创建第二个三角形，或不超出范围，并确保你的主函数工作正常。
4. 创建一个L型域的图像（在step-1中添加一个函数`third_grid`），并使用`GridOut::write_vt()'函数将其输出为`third_grid.vtk'。
5. 围绕再入角自适应地细化L型网格三次（在你已经做过的全局细化之后），但有一个变化：细化所有单元，其中心与再入角之间的距离小于 小区中心与再入角之间的距离小于 2/3。
6. 创建一个辅助函数，接收一个三角形的引用，并返回一个元组，其内容如下：
   返回一个包含以下信息的元组：层数、单元数、活动单元数。用你所有的网格测试这个函数。
7. 取消对`triangulation.reset_manifold(0)`行的注释。
   请问，`second_grid()`一行。现在会发生什么？把你的发现作为注释写在文件中的那一行下面，在你做完实验后，把命令注释掉。
8. 额外任务：创建一个代表圆环表面的网格，并在全局范围内全局性地对其进行2次细化。输出到vtk格式并检查输出结果。注意你的三角形需要是 "Triangulation<2,3>"的类型。
9. 额外任务：看一下step-49，在你修改过的step-1程序中阅读包含的.msh文件。

### Lab-02: step-2

 1. 阅读step-2的文档： <https://www.dealii.org/current/doxygen/deal.II/step_2.html>。
 
 2. 运行step-2。在浏览器中看一下稀疏模式。
 3. 如果你把多项式的度数从1增加到2或3，模式会有什么变化？
 4. 如果你使用全局细化的（比如3倍）单位方，模式会有什么变化？
 5. 这些模式是对称的吗？为什么/为什么不是？
 6. 对于一个Q1元素，你期望每行有多少个条目？Q1元素（假设每个顶点周围有四个单元）？检查一下 对b)中的网格来说是真实的(寻找`row_length(i)`并输出每一行的数据)。每一行的长度）。) 你能不能构建一个二维网格（没有悬挂的 节点），并且有更多条目的行？
 7. 假设每个顶点周围有四个单元，那么在稀疏模式中每行有多少个条目是关于Q2和 假设每个顶点周围有四个单元，那么Q2和Q3元素每行有多少条目？
 8. 打印第42行的所有条目，即原始重新编号的稀疏模式。模式的所有条目。
 9. 奖励：计算并输出统计数据，如未知数的数量、稀疏模式的带宽、平均数等。未知数的数量，稀疏模式的带宽，每行的平均条目数和填充率。每行的平均条目数，以及填充率。


# Lab 03 - Poisson问题

## 一般说明

对于下面的每一点，用函数来扩展`step-3`类，执行指定的任务，尽量减少你复制的代码量 尽量减少复制和粘贴的代码量，可能的话，通过向现有函数添加参数来重组现有代码 函数，并生成类似于`run`方法的封装器（例如`run_exercise_3`)

一旦你创建了一个执行给定任务的函数，将其添加到 `step-3-test.cc`文件，并确保所有的练习都通过 `gtest`可执行文件，例如，为每个练习添加一个测试，如 下面的片段。

```C++
TEST_F(Step3Tester, Exercise3) {
   run_exercise_3();
}
```
在本实验室结束时，你将有一个代码，可以解决二维的泊松问题，在不同的领域类型上，有任意的右手边，任意的狄里奇边界条件，可变的有限元程度，以及不同的细化水平。

该程序将为其执行读取参数文件。这是简单但可扩展的有限元应用的最小起点，并且已经触及了一般和可扩展的有限元代码的许多重要和基本特征。 可扩展的有限元代码的许多重要的基本特征。

## Lab-03 

### step-3

1.  参见步骤3的文档 <https://www.dealii.org/current/doxygen/deal.II/step_3.html>

2.  编译并运行step-3。检查源文件和头文件。将所有`#include<...>`指令，你可以从头文件`step-3.h`移到源文件`step-3.cc`。

3.  打开vtk输出，在paraview或visit中查看解决方案。弄清楚如何通过解决方案的变量利用“warp”查看方案。保存一张图片在 "figures "子目录下保存你的可视化图片，命名为 "step-3-0.png"。并提交到你的存储库。

4.  按照教程描述中 "修改边界条件的类型 "的指示中的说明，绘制你的输出，并将图保存在绘制输出结果，并在 "figures "子目录下保存一个图，命名 "step-3-1.png"，并提交到你的资源库中。储存库。确保测试是通过`gtest`运行的。

5.  现在也做 "最后一点的轻微变化"，但使用数值-0.5作为指标1的边界。将你的输出保存为`step-3-2.png`。并将其提交到版本库的`figures`目录下。请确保该测试是通过`gtest`运行的。

6.  改变设置，使$f=0$。将你的输出保存为`step-3-3.png`。并提交到版本库中的`figures`目录。请确保该测试是通过`gtest`运行的。

7.  7.切换到一个L型域，并尝试将迪里希特和诺伊曼的组合。Dirichlet和Neumann边界条件。通过实验，确定识别与再入角相邻的面，并只在那里应用迪里希特条件。只适用于那里。将你的输出保存为`step-3-4.png`。并将其提交到资源库中的`figures`目录。请确保该测试是通过`gtest`运行的。

8.  奖励：做 "均值收敛"（通过阅读的文档）。 你能看到阶数 $h^2$ 吗？增加多项式的阶数（你需要增加程序中所有阶数的你需要增加程序中的所有阶数！）并检查平均数的收敛的平均值，并确保测试是通过`gtest`运行的。

### ParameterHandler

1.  按照`ParameterAcceptor'类的文件规定  (来自 "这个类的典型用法")。从`ParameterAcceptor`派生出你的`Step3`类。从`ParameterAcceptor`派生出你的`Step3`类，并使用`ParameterAcceptor`类中的字符串`"Step3 "初始化。字符串`"Step3"`。添加一个函数来初始化这个类，给定一个参数文件名，即`initialize("parameters.prm")'，通过调用 `ParameterAcceptor::initialize(...)`。
    
2.  在你的 "Step3 "类中添加以下参数，并确保它们在整个代码中被使用。 在你的代码中使用
      - 有限元程度
      - 全局细化的数量
      - 输出文件名

    请注意，`FE_Q`类需要在*初始化有限元度变量之后建立。初始化有限元素度变量之后。你可以在使用前创建你可以在使用它之前创建它，然后通过以下方式引用它 `DoFHandler::get_fe()`，然后在以后需要时引用它。

3.  在你的 "Step3 "类中添加两个 "FunctionParser<2>"成员，一个用于边界条件，另一个用于强制项。添加两个`std::string`成员，代表他们的数学表达式，和一个`std::map<std::string, double>`代表在数学表达式中使用的常数，并确保它们被添加为你的类的参数，作为
      - 强制项表达式
      - 边界条件表达式
      - 问题常数

    并用给定的表达式和给定的常数初始化这两个`FunctionParser<2>`对象。确保你的`Step3'类使用这两个函数，在那里计算边界条件，并在那里组装右手边。

4.  按照`步骤-70'的思路，在`Step3'类中添加两个额外的参数，以便在运行时可以选择边界条件。
    `Step3`类中增加两个参数，以便在运行时选择创建什么样的网格，以及用什么样的参数。作为参数使用
      - 网格生成器功能
      - 网格生成器参数

5.  用`parameters`子目录下的所有文件测试你的应用程序。(通过`gtest`应用程序运行测试)，并为每个给定的参数创建一个 为每个给定的参数创建一个可视化的文件，其名称与 参数文件的名称，也就是说，运行测试的可视化文件是 `parameters/hyper_shell.prm`应该保存在一个名为`figures/hyper_shell.png`的图片上。提交所有你生成的数字。

    注意，对于老版本的google测试，你可能需要单独运行每个 测试（这是google testsuite和`ParameterAcceptor`之间的一个已知问题）。`ParameterAcceptor`类之间的问题）。) 要运行一个单独的测试，你可以使用 `--gtest_filter`命令行选项，即以下命令

    ```
    ./gtest --gtest_filter=Step3Tester.HyperShell
    ```

    将只运行测试`Step3Tester.HyperShell'。

    如果你想一个接一个地运行所有测试，你可以调用`ctest`工具。实用程序，它将执行`gtest`命令，通过gtest过滤器对你提供的每个测试进行测试。的测试，例如，产生类似于以下的输出 的输出。

    ```
    dealii@66c3be678611:/workspace/sissa-mhpc-lab-03/build-container$ ctest
    测试项目 /workspace/sissa-mhpc-lab-03/build-container
        开始1：Step3Tester.MakeGrid
    1/6 测试#1：Step3Tester.MakeGrid .............   通过 0.68 秒
        开始 2: Step3Tester.HyperCube
    2/6 测试 #2: Step3Tester.HyperCube ............   通过0.78秒
        开始3：Step3Tester.HyperShell
    3/6 测试 #3: Step3Tester.HyperShell ...........   通过 1.90 秒
        开始4：Step3Tester.HyperBall
    4/6 测试 #4: Step3Tester.HyperBall ............   通过 1.28 秒
        开始 5: Step3Tester.SinRhs
    5/6 测试 #5: Step3Tester.SinRhs ...............   通过 0.79 秒
        开始 6: Step3Tester.SinBc
    6/6 测试 #6: Step3Tester.SinBc ................   通过 1.26 秒
    
    100%的测试通过，6个测试中有0个测试失败
    ```

# Lab 04 - 模板参数和收敛率


## 说明

对于下面的每一点，用执行指定任务的函数来扩展`Poisson`类，尽量减少你复制和粘贴的代码量，可能的话通过向现有函数添加参数来重组现有代码，并生成类似于`run`方法的包装器（如 `run_exercise_3`）。

一旦你创建了一个执行给定任务的函数，将其添加到`poisson-tester.cc`文件中，并确保所有练习都通过`gtest`可执行程序运行，例如，为每个练习添加一个测试，如下所示片段。

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

在本实验室结束时，你将有一个代码来解决任意维度的泊松问题。泊松问题的代码，该代码可以在不同的领域类型、不同的边界条件下，用任意程度的拉格朗日有限元解决任意维度的泊松问题。度的拉格朗日有限元，在不同的领域类型上，用不同的边界条件和 不同的函数来定义右手边、刚度系数和强制项。僵化系数和强迫项的定义。

该问题将在连续的细化网格上运行，我们将验证 不同阶数的拉格朗日有限元空间的Bramble-Hilbert法则。使用Python建立人造解决方案，并使用Latex、Tikz 使用latex、tikz和pgfplots绘制误差收敛表。

该程序将建立在你对Step-3的实现之上，借鉴 step-4"、"step-5 "和 "step-7"。

## Lab-04 

### step-4

1.  参见step-4的文档 <https://www.dealii.org/current/doxygen/deal.II/step_4.html>

2.  编译并运行step-4。检查源文件和头文件。

3. 将`lab-03'中`步骤3'的实现复制到文件中 `source/poisson.cc`, `include/poisson.h`, 和 `tests/poisson-tester.cc`, 确保你的所有文件和类都正确地重命名。 确保你将所有的文件和类正确地重命名为 "Poisson"。

4. 在你的 "Poisson "类中添加模板参数"<int dim>"，以`step-4'为例，确保你的程序在2D和3D中都能正确运行。确保你的程序在2D和3D中都能正确运行。

5. 添加参数
   
    - 细化周期的数量 `Number of refinement cycles`
    - 精确解表达式 `Exact solution expression`
   
   和相应的成员变量（即`n_cycles`, 确切的解决方案表达式 "和"确切的解决方案"），并对每个细化周期再次运行泊松问题。对每个细化周期运行一次全局细化的泊松问题，确保输出每个细化的结果。确保你以`vtu`格式分别输出每个细化周期的结果。格式输出，即如果 "输出文件名 "是 "poisson_2d"，"细化周期数"是3，那么细化周期数是3，你应该输出

     - `poisson_2d_0.vtu`
     - `poisson_2d_1.vtu`
     - `poisson_2d_2.vtu`

    其中`poisson_2d_0.vtu`中的解决方案应该有`全局细化数`细化，`poisson_2d_1.vtu`应该有`全局细化数`+1细化，而`poisson_2d_2.vtu`应该有`全局细化数`+2细化。

6. 在你的 "泊松 "类中添加一个 "ParsedConvergenceTable "对象（见https://www.dealii.org/current/doxygen/deal.II/classParsedConvergenceTable.html)
并在参数文件的 "错误表 "小节中添加其参数。即在`Poisson`构造函数中添加以下几行代码。
```
this->prm.enter_subsection("Error table");
error_table.add_parameters(this->prm);
this->prm.leave_subsection();
```

7. 设置边界条件、强制函数和精确表达式以得到
u(x,y)=sin(pi*x)*cos(pi*y)"。添加一个方法`compute_error()`到`Poisson`类中，调用`ParsedConvergenceTable::error_from_exact`方法和你上面创建的`exact_solution`函数。请确保你将L2和H1的错误以文本格式输出到一个文件中。文本格式输出到一个文件。使用`jupyter'笔记本
`manufactured_solutions.ipynb`来构建非微妙的精确解。

8. 添加参数 "刚度系数表达式 "和相应的的成员，这样你要解决的问题就是 $-div(coefficient(x)\nabla u) = f(x)$。

9. 构建一个(非星形的!)人造解，其中`系数`是一个不连续的函数。不连续的函数。请注意，制造的解决方案可能需要一个注意，制造的解决方案可能需要在右边有一个不连续的强制项，但不应该有其他的
的奇异性，也就是说，你需要确保$coefficient(x)\nabla u$是连续的，也就是说，$\nabla u$有一个跳跃取决于 "系数 "的跳变。输出误差表，并对你观察到的错误率进行评论。当增加多项式的阶数时，情况是否有所改善？多项式阶数时，情况是否有所改善？为什么？

10. (可选)使用latex子目录中提供的latex文件来为你的解决方案生成专业的收敛图。


# Lab 05 - 边界条件和约束

### 说明

对于下面的每一点，用函数扩展 "Poisson "类，以完成指定的任务。
执行指定的任务，尽量减少你复制和粘贴的代码量。尽量减少复制和粘贴的代码量，可能的话，通过给现有的函数添加参数来重组现有的代码。并生成类似于 "运行 "方法的封装器（如`run_exercise_3`）

一旦你创建了一个执行给定任务的函数，将其添加到`poisson-tester.cc`文件，并确保所有练习都通过
`gtest`可执行程序运行，例如，为每个练习添加一个测试，如下所示片段。

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

在本实验结束时，你将修改你的Poisson代码以便允许在域的不同部分有非均质的诺伊曼边界条件。
你将为求解器添加更多的选项，使之能够使用直接求解器，或使用其他方法。使用直接求解器，或一些更复杂的预处理器。

## Lab-05

### step-5

1. 添加参数
   
    - `Neumann boundary condition expression`
    - `Dirichlet boundary ids`
    - `Neumann boundary ids`
   
   和相应的成员变量到你的泊松问题问题上，使用 `std::set<dealii::types::boundary_id>`为最后两个参数。

2. 修改你的Dirichlet边界条件的实现，以便将Dirichlet函数应用于参数文件中指定的所有边界ID。
文件中指出的所有边界ID。

3. 使用上面定义的函数实现诺伊曼边界条件。在参数文件中指出的纽曼边界的id上实现纽曼边界条件。

4. 创建一个函数，从一个字符串中解析参数，在测试基础设施中使用。创建一个函数，从字符串中解析参数，用于测试基础设施

5. 创建一些测试，在非常小的网格上实际解决泊松问题。用非常简单但非琐碎的边界条件组合，在不同的域上实际解决泊松问题。不同的领域，并验证你的发现的正确性。例如
包括使用全局线性（二次）精确函数，用线性二次方）有限元，并验证你所做的误差实际上是零（在这种情况下，全局插值给出了精确解，因此有限元也应该提供精确的解决方案)

6. 修改你的代码，使用 "AffineConstraints "代替`VectorTools::apply_boundary_values`（见`步骤6`的文档）。再次运行所有的测试，并验证你的新代码是否通过了自己的检查。

7. 添加参数。
    - 本地预精化网格尺寸表达式
   
   和相应的成员。在创建网格时，不是简单地全局细化的固定次数（由参数 
   全局细化次数"），而是在局部细化你的网格，当函数在一个单元格的中心被评估为大于实际的单元格直径。请确保在局部细化的次数等于你所做的局部细化循环数等于参数`全局细化次数`。将上述函数设置为 "0 "应产生与之前相同的结果，即固定的全局细化次数，直至 "全局细化次数"。

8.  确保你正确计算了悬挂节点的约束，并且你的解算器也能在悬空节点上工作。

9. 给你的求解器添加一个预处理程序。如果你安装了trilinos并且它被配置在`deal.II`中，使用它的代数多网格预处理程序（它也适用于`deal.II`）。它也适用于`deal.II`的矩阵）。 否则就使用库中其他可用的先决条件，否则使用库中其他可用的先决条件。验证你的求解器现在是否更快。


# Lab 06 - 后验误差估计和自适应有限元分析

### 说明

对于下面的每一点，用函数扩展 "Poisson "类，以便执行指定的任务，尽量减少你复制和粘贴的代码量。
尽量减少复制和粘贴的代码量，可能的话，通过给现有的函数添加参数来重组现有的代码。并生成类似于 "运行 "方法的封装器（如`run_exercise_3`）。

一旦你创建了一个执行给定任务的函数，将其添加到`poisson-tester.cc`文件，并确保所有练习都通过
`gtest`可执行程序运行，例如，为每个练习添加一个测试，如下所示片段。

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

在本实验结束时，你将修改你的Poisson代码以便允许在域的不同部分有非均质的诺伊曼边界条件。
你将为求解器添加更多的选项，使之能够使用直接求解器，或使用其他方法。使用直接求解器，或一些更复杂的预处理器。

## Lab-06
### step-6

1.  参见步骤6的文档 <https://www.dealii.org/current/doxygen/deal.II/step_6.html>

2. 添加参数
   
    - `Mapping degree`
    - `Marking strategy`
    - `Estimator type`
    - `Coarsening and refinement factors`
   
    其中`映射程度`控制代码中使用的映射程度。标记策略 "是在 "全局|固定分数|固定数字 "之间选择。
估算器类型 "可在 "精确 "或 "凯利 "或 "残留 "之间选择，以及
粗化和细化因子 "是一个 "std::pair<double, double>"，包含传递给 "std::down "的参数。
参数传递给`GridRefinement::refine_and_coarsen_fixed_*`函数。

3. 确保所有的 "FEValues "类都使用正确顺序的映射，并且确保你在输出中也使用正确的映射。

4. 在`Poisson`中添加一个`Vector<float`字段`error_per_cell`，由`estimate`方法填充。

5. 增加一个方法`residual_error_estimator`，使用`FEInterfaceValues`和`FEValues`计算残差估计值。

6. 在 "Poisson "类中增加一个方法 "estimate"，用于计算准确值和计算值之差的H1 seminorm。
如果 "估算器类型 "为 "精确"，则计算精确解和计算解之间的差异的H1 seminorm。是 "exact"，如果 "Estimator类型 "是 "kellyErrorEstimator<dim>::impair"，则调用"KellyErrorEstimator<dim>::impair"。如果 "估计器类型 "是 "kelly"，则调用 "residual_error_estimator"，如果 "估计器类型 "是 "residual"，则调用 "residual_error_estimator"。

7. 在收敛表中也加入你计算的估计器。这应该是与`Estimator type`是`exact`的情况下的`H1`半正态相同。

8. 以L型域上计算出的精确解为例，计算自适应有限元的速率计算自适应有限元方法在自由度数方面的收敛率。使用上述三个估计器的自由度收敛率

9. 设定边界条件为零，强迫项等于4，精确解等于`-x^2-y^2+1`，并在圆心在原点、半径为1的圆上求解问题。在一个中心在原点、半径为1的圆上，对各种有限元和映射度进行求解。当映射度不高时，你会观察到什么？
当映射度与有限元度不一致时，你观察到什么？你如何解释这个问题？

10. 创建一个测试，完全重现 "第6步 "的行为，仅使用你的参数文件


# Lab 07 - Adaptive FEM and shared memory parallelization

### 说明

对于下面的每一点，用函数来扩展"Poisson"类，以便执行指定的任务，尽量减少你复制和粘贴的代码量。尽量减少复制和粘贴的代码量，可能的话，通过给现有的函数添加参数来重组现有的代码。函数，并生成类似于 "运行 "方法的封装器（如`run_exercise_3`）。

一旦你创建了一个执行给定任务的函数，将其添加到`poisson-tester.cc`文件，并确保所有练习都通过
`gtest`可执行程序运行，例如，为每个练习添加一个测试，如下所示片段。

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

在本实验结束时，你将修改你的Poisson代码，以便在多线程上使用共享内存并行化来运行在多线程上使用共享内存并行化运行，并且对基于任务的并行化有一定了解。你将对基于任务的并行化有一些了解。

### Step-7

1. 用一个 "TimerOutput "类装备你的 "泊松 "求解器，并提取 " Poisson "类的每个方法的时间信息（见https://www.dealii.org/current/doxygen/deal.II/classTimerOutput.html）使用范围内的定时器。进行几次测试运行，并把你用来运行代码的参数你用来运行代码的参数文件，以及计时结果。

2. 阅读关于共享内存并行化的文件。
   https://www.dealii.org/current/doxygen/deal.II/group__threads.html

3. 添加参数
   
    - 最大的线程数
   
    允许使用一个整数。值为`-1`意味着*选择自动*。对于任何其他的数字，确保你调用`MultithreadInfo::set_thread_limit`，并给出参数。

4. 使用`Threads::Thread`或`Threads::Thread`进行*粗粒度*并行化。在你的 "泊松 "求解器的每个主要任务上使用 "Threads::ThreadGroup "进行粗粒度*并行化。实验改变上面的参数，并报告你的代码在使用多线程时的加速情况。并报告使用多线程时的代码速度，与上面第1步中得到的结果相比。

5. 使用`Threads::split_range`将系统与尽可能多的CPUS并行组装。你有多少CPUS就有多少。检查你是否确实得到了改善。

6. 使用deal.II默认的`ScratchData`和`CopyData`对象（即。
https://www.dealii.org/current/doxygen/deal.II/classMeshWorker_1_1ScratchData.html
和https://www.dealii.org/current/doxygen/deal.II/structMeshWorker_1_1CopyData.html
将`Threads::split_range`组件替换为基于使用`WorkStream::run`方法，并将运行时间与你原来的代码进行比较。

7. 使用`Doxygen`语法记录`Poisson`类的所有方法和成员。确保GitHub行动`documentation.yaml`是有效的，并且自动生成你的代码的文档。如果一切正常。你应该能够以`https://username.github.io/sissa-mhpc-lab-07-username/`，一旦动作完成后，你就可以以‘web’的形式访问它。