# deal.ii-course-practice

# Theory and Practice of Finite Element Methods

## Laboratory \#1

Goal of the laboratory: 
- set up a modern programming environment
- familiarse with git, VSCode, CMake, GoogleTest, and C++

Tasks:
- split the source file `pythagoras.cc` into two files, one containing
the `main` function, and the other one containing the tests
- move them to a `source` directory, and name them `main.cc` and `point_tests.cc`
- modify the `CMakeLists.txt` file to ensure the program still compiles and run
- Add three more tests that checks that the three-dimensional version of the `Point<dim>` class works as expected


#  Lab 02 - Triangulation, DoFHandler, FiniteElement
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

###  Setup

- You will need the following programs:
  - Visual Studio Code with the `Remote Development` plugin installed
  - Paraview (ver > 5.4 would be better)
-  Check that the `.devcontainer` directory contains a configuration that works foryou. When you open the root directory of the repository, it should prompt you toreopen in container. If not, try explictily with the command `F1 -> Reopen folder incontainer`
- If opening in container worked, you should be able to run a terminal (`CTRL + ~`)and see something like
```
To run a command as administrator (user "root"), use "sudo <command>".
See "man sudo_root" for details.
dealii@3f25621373bd:/workspaces/02_tria_dh_fem$ 
```
- The folder contains a `CMakeLists.txt` file. VSCode should detect this, and providecommands to configure and build your application, at the bottom of the terminal window


### Setup if you have ubuntu/debian and DO NOT want to use Docker

In this case you can copy the commands from the docker image that you can find here:

- <https://github.com/dealii/docker-files/blob/master/dependencies-focal/Dockerfile>

to install all the dependencies. After this, you can install the latest version 
of the library using the commands you find here:

- <https://github.com/dealii/dealii/blob/master/contrib/docker/Dockerfile>

### Lab-02: step-1

1.  Read the documentation of step-1 at
    -   <https://www.dealii.org/current/doxygen/deal.II/step_1.html>
2. Compile and run `step-1` inside the VSCode container and look at the output.
3. Understand what happens if you now call `first_grid(triangulation)` and `second_grid(triangulation)` one after the other. Fix the problem you get, **without** creating a second triangulation, or without going out of scope, and make sure your main function works correctly.
4. Create an image of an L-shape domain (add a function `third_grid` to step-1) 
   with one global refinement, and output it as `third_grid.vtk` using the `GridOut::write_vt()` function.
5. Refine the L-shaped mesh adaptively around the re-entrant corner
   three times (after the global refinement you already did), but with a twist: refine all cells with the distance
   between the center of the cell and re-entrant corner is smaller than
   2/3.
6. Create a helper function that takes a reference to a Triangulation and 
   returns a tuple with  the following information: number of levels, number of cells, number of active cells. Test this with all of your meshes.
7. Un-comment the `triangulation.reset_manifold(0)` line in
   `second_grid()`. What happens now? Write your findings as a comment under the line in the file, and leave the command commented out after you experimented with it.
8. Bonus: Create a mesh that represents the surface of a torus and refine
   it 2 times globally. Output to vtk format and check the output. Note
   that your Triangulation needs to be of type ``Triangulation<2,3>``,
   which we will discuss later this week.
9. Bonus: Take a look at step-49 and read the included .msh file in your modifiedstep-1 program.

### Lab-02: step-2

 1. Read the documentation of step-2 at
    -  <https://www.dealii.org/current/doxygen/deal.II/step_2.html>

 1. Run step-2. Look at the sparsity patterns in your browser.

 2. How does the pattern change if you increase the polynomial degree from 1 to 2 orto 3?

 3. How does the pattern change if you use a globally refined (say 3 times) unitsquare?

 4. Are these patterns symmetric? Why/why not?

 5. How many entries per row in the sparsity pattern do you expect for a
    Q1 element (assuming four cells are around each vertex)? Check that
    this is true for the mesh in b) (look for `row_length(i)` and output
    them for each row). Can you construct a 2d mesh (without hanging
    nodes) that has a row with more entries?

 6. How many entries per row in the sparsity pattern are there for Q2 and
    Q3 elements, again assuming four cells around each vertex?

 7. Print all entries for row 42 for the original renumbered sparsity
     pattern.

 8. Bonus: Compute and output statistics like the number of
    unknowns, bandwidth of the sparsity pattern, average number of
    entries per row, and fill ratio.


#  Lab 03 - Poisson Problem
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

#  Lab 03 - Poisson Problem
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

## General Instructions

For each of the point below, extend the `step-3` class with functions that 
perform the indicated tasks, trying to minimize the amount of code you copy
and paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g., 
`run_exercise_3`).

Once you created a function that performs the given task, add it to the 
`step-3-tester.cc` file, and make sure all the exercises are run through
the `gtest` executable, e.g., adding a test for each exercise, as in the 
following snippet: 

```C++
TEST_F(Step3Tester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have a code that solves a Poisson 
problem in two dimensions, on different domain types, with arbitrary right hand
side, arbitrary Dirichlet boundary conditions, variable finite element degree, 
and different levels of refinement.

The program will read parameter files for its execution. This is the minimal 
starting point for simple but extensible Finite Element applications, and 
touches already many important and fundamental characteristics of general and 
extensible Finite Element codes.

## Lab-03 

### step-3

1.  See documentation of step-3 at
    <https://www.dealii.org/current/doxygen/deal.II/step_3.html>

2.  Compile and run step-3. Examine the source and header files. Move all 
    `#include<...>` directives that you can from the header `step-3.h` to
    the source file `step-3.cc`.

3.  Open the vtk output and visualize the solution in paraview or visit. 
    Figure out how to warp the solution by the solution variable. Save a picture
    of your visualization on the `figures` subdirectory, named `step-3-0.png`, 
    and commit it to your repository.

4.  Follow the instructions in “Modify the type of boundary condition”
    in the description of the tutorial, plot your output, and save a figure in
    the `figures`  subdirectory, named `step-3-1.png`, and commit it to your
    repository. Make sure the test is run through `gtest`.

5.  Now also do “A slight variation of the last point” but use the value
    -0.5 for the boundary with indicator 1. Save your output as `step-3-2.png`, 
    and commit it to the `figures` directory in your repository. Make sure the
    test is run through `gtest`.

6.  Change the setup to have $f=0$. Save your output as `step-3-3.png`, 
    and commit it to the `figures` directory in your repository. Make sure the
    test is run through `gtest`.

7.  Switch to an L-shaped domain and experiment with a combination of
    Dirichlet and Neumann boundary conditions. By experimentation, identify
    the faces adjacent to the re-entrant corner and apply Dirichlet conditions
    only there. Save your output as `step-3-4.png`, 
    and commit it to the `figures` directory in your repository. Make sure the
    test is run through `gtest`.

8.  Bonus: Do “Convergence of the mean” (read through the end of the 
    documentation of ). Can you see the order $h^2$?
    Increase the polynomial order (you need to increase all orders of
    the quadratures in the program!) and check the convergence of the
    mean now. Make sure the test is run through `gtest`.

### ParameterHandler

1.  Follow the documentation of the class `ParameterAcceptor` 
    (from "A typical usage of this class"). Derive your `Step3` class from 
    `ParameterAcceptor`, and initialize the `ParameterAcceptor` class using 
    the string `"Step3"`. Add a function that initializes the class given a 
    parameter file name, i.e., `initialize("parameters.prm")`, by calling 
    `ParameterAcceptor::initialize(...)`
    
2.  Add the following parameters to your `Step3` class, and make sure they 
    are used throughout your code

      - Finite element degree
      - Number of global refinements
      - Output filename

    Notice that the `FE_Q` class will need to be built *after* the 
    initialization of the Finite element degree variable. You can create
    it right before you use it, and then reference to it by 
    `DoFHandler::get_fe()` when you need it later.

3.  Add two `FunctionParser<2>` members to your `Step3` class, one for the 
    boundary conditions, and one for the forcing term. Add two `std::string` 
    members, representing their mathematical expressions, and a 
    `std::map<std::string, double>` representing the constants to use in your
    mathematical expression, and make sure they are added as parameters of your 
    class as
    
      - Forcing term expression
      - Boundary condition expression
      - Problem constants
    
    and initialize the two `FunctionParser<2>` objects with the given 
    expressions and the given constants. Make sure your `Step3` class uses 
    these two functions where the boundary conditions are computed and where
    the right hand side is assembled

4.  Following the ideas in `step-70`, add two additional parameters to the 
    `Step3` class to select at run time what grid to create, and with what 
    arguments. Use as parameters

      - Grid generator function
      - Grid generator arguments

5.  Test your application with all the files in the `parameters` subdirectory,
    (running the test through the `gtest` application), and create a 
    visualization for each of the given parameters, with the same
    name of the parameter file, i.e., the visualization of the test that runs
    `parameters/hyper_shell.prm` should be saved on an image named `figures/hyper_shell.png`. Commits all your generated figures.
    Notice that for older versions of google test, you may need to run each
    test individually (this is a known issue between google testsuite and the
    `ParameterAcceptor` class). To run an individual test, you can use the 
    `--gtest_filter` command line option, i.e., the following command

    ```
    ./gtest --gtest_filter=Step3Tester.HyperShell
    ```

    will run only the test `Step3Tester.HyperShell`.

    If you want to run all tests one after the other you can call the `ctest` 
    utility, which will execute the `gtest` command passing the gtest filter for
    each of the tests you provided, producing, for example, an output similar 
    to:

    ```
    dealii@66c3be678611:/workspace/sissa-mhpc-lab-03/build-container$ ctest
    Test project /workspace/sissa-mhpc-lab-03/build-container
        Start 1: Step3Tester.MakeGrid
    1/6 Test #1: Step3Tester.MakeGrid .............   Passed    0.68 sec
        Start 2: Step3Tester.HyperCube
    2/6 Test #2: Step3Tester.HyperCube ............   Passed    0.78 sec
        Start 3: Step3Tester.HyperShell
    3/6 Test #3: Step3Tester.HyperShell ...........   Passed    1.90 sec
        Start 4: Step3Tester.HyperBall
    4/6 Test #4: Step3Tester.HyperBall ............   Passed    1.28 sec
        Start 5: Step3Tester.SinRhs
    5/6 Test #5: Step3Tester.SinRhs ...............   Passed    0.79 sec
        Start 6: Step3Tester.SinBc
    6/6 Test #6: Step3Tester.SinBc ................   Passed    1.26 sec
    
    100% tests passed, 0 tests failed out of 6
    ```


    #  Lab 04 - Template parameters and convergence rates
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

## General Instructions

For each of the point below, extend the `Poisson` class with functions that
perform the indicated tasks, trying to minimize the amount of code you copy and
paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g.,
`run_exercise_3`).

Once you created a function that performs the given task, add it to the
`poisson-tester.cc` file, and make sure all the exercises are run through the
`gtest` executable, e.g., adding a test for each exercise, as in the following
snippet:

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have a code that solves a Poisson
problem in arbitrary dimensions, with Lagrangian finite elements of arbitrary
degree, on different domain types, with different boundary conditions, and
different functions for the definition of the right hand side, the stiffness
coefficient, and the forcing term.

The problem will run on successively refined grids, and we will verify
Bramble-Hilbert lemma for Lagrangian finite element spaces of different order,
building manufactured solutions using python, and plotting error convergence
tables using latex, tikz, and pgfplots.

The program will build on top of your implementation of Step3, drawing from
`step-4`, `step-5`, and `step-7`.

## Lab-04 

### step-4

1.  See documentation of step-4 at
    <https://www.dealii.org/current/doxygen/deal.II/step_4.html>

2.  Compile and run step-4. Examine the source and header files.

3. Copy your implementation of `step-3` from `lab-03` to the files
`source/poisson.cc`, `include/poisson.h`, and `tests/poisson-tester.cc`, make
sure you rename correctly all your files and classes to `Poisson`.

4. Add the template parameter `<int dim>` to your `Poisson` class, following
`step-4` as an example, and make sure that your program runs correctly both in
2D and in 3D.

5. Add the parameters
   
    - `Number of refinement cycles`
    - `Exact solution expression`
   
   and the corresponding member variables (i.e., `n_cycles`,
   `exact_solution_expression`, and `exact_solution`) and run the Poisson
   problem again for each refinement cycle with one global refinement, making
   sure you output the result for each refinement cycle separately in `vtu`
   format, i.e., if `Output filename` is `poisson_2d`, and `Number of
   refinement cycles` is 3, you should output

     - `poisson_2d_0.vtu`
     - `poisson_2d_1.vtu`
     - `poisson_2d_2.vtu`

    where the solution in `poisson_2d_0.vtu` should have `Number of global
    refinements` refinements, `poisson_2d_1.vtu` should have `Number of global
    refinements` +1 refinements, and `poisson_2d_2.vtu` should have `Number of
    global refinements` +2 refinements.

6. Add a `ParsedConvergenceTable` object to your `Poisson` class (see
https://www.dealii.org/current/doxygen/deal.II/classParsedConvergenceTable.html)
and add its parameters in the subsection `Error table` of the parameter file,
i.e., in the `Poisson` constructor add the following lines of code:
```
this->prm.enter_subsection("Error table");
error_table.add_parameters(this->prm);
this->prm.leave_subsection();
```

7. Set the boundary conditions, forcing function, and exact expression to get
the manufactured solution `u(x,y)=sin(pi*x)*cos(pi*y)`. Add a method
`compute_error()` to the `Poisson` class, that calls the
`ParsedConvergenceTable::error_from_exact` method with the `exact_solution`
function you created above. Make sure you output both the L2 and H1 error in
text format to a file. Play with the `jupyter` notebook
`manufactured_solutions.ipynb` to construct non-trivial exact solutions.

8. Add a parameter `Stiffness coefficient expression` and the corresponding
members to the `Poisson` class, so that the problem you will be solving is 
$-div(coefficient(x)\nabla u) = f(x)$.

9. Construct a (non-singular!) manufactured solution where `coefficient` is a
discontinuous function. Notice that the manufactured solution may need a
discontinuous forcing term on the right hand side, but should not have other
types of singularities, that is, you need to make sure that
$coefficient(x)\nabla u$ is continuous, i.e., that $\nabla u$ has a jump
depending on the jump of the `coefficient`. Output the error tables, and
comment on the error rates you observe. Do things improve when increasing the
polynomial order? Why?

10. (optional) Use the latex file provided in the latex subdirectory to
generate professional convergence plots for your solutions.


#  Lab 05 - Boundary conditions and constraints
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

## General Instructions

For each of the point below, extend the `Poisson` class with functions that
perform the indicated tasks, trying to minimize the amount of code you copy and
paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g.,
`run_exercise_3`).

Once you created a function that performs the given task, add it to the
`poisson-tester.cc` file, and make sure all the exercises are run through the
`gtest` executable, e.g., adding a test for each exercise, as in the following
snippet:

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have modified your Poisson code to
allow also non-homogeneous Neumann boundary conditions on different parts of
the domain, and you will have added some more options to the solver, enabling
usage of a direct solver, or of some more sofisticated preconditioners.

## Lab-05

### step-5

1. Add the parameters
   
    - `Neumann boundary condition expression`
    - `Dirichlet boundary ids`
    - `Neumann boundary ids`
   
   and the corresponding member variables to your Poisson problem problem, using 
   `std::set<dealii::types::boundary_id>` for the last two parameters

2. Modify your implementation of Dirichlet boundary conditions, in order to
apply the Dirichlet function to all boundary ids indicated in the parameter
file

3. Implement Neumann boundary conditions, using the function defined above, on
the ids of the Neumann boundary indicated in the parameter file

4. Create a function that parses parameters from a string, to be used in the
testing infrastructure

5. Create a few tests that actually solve a Poisson problem on very small grids
with very simple but non-trivial combinations of boundary conditions, on
different domains, and verify the correctness of your findings. Examples
include using globally linear (quadratic) exact functions, with linear
(quadratic) finite elements, and verify that the error you make is actually
zero (in this case, global interpolation gives the exact solution, therefore
the finite element should also provide the exact solution)

6. Modify your code to use `AffineConstraints` instead of
`VectorTools::apply_boundary_values` (see the documentation of `step-6`). Run
again all tests, and verify that you pass your own checks again with the new
code

7. Add the parameters:
    - `Local pre-refinement grid size expression`
   
   and the corresponding members. When creating the grid, instead of simply
   refining globally a fixed number of times (given by the parameter 
   `Number of global refinements`), refine locally your grid when the function
   above evaluated in the center of a cell is larger then the actual cell 
   diameter. Make sure you stop refining locally if the number of local
   refinement cycles you did is equal to the parameter 
   `Number of global refinements`. Setting the function above to `0` should produce the same results as before, i.e., a fixed number of global refinements up to `Number of global refinements`.

8. Make sure you compute correctly the hanging node constraints, and that your
solver works also with hanging nodes

9. Add a preconditioner to your solver. If you have trilinos installed and it
is configured inside `deal.II`, use its algebraic multigrid preconditioner (it
works also with `deal.II` matrices). Otherwise use one of the other available
preconditioners in the library. Verify that your solver is now faster.


#  Lab 06 - A posteriori error estimation and adaptive FEM
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

## General Instructions

For each of the point below, extend the `Poisson` class with functions that
perform the indicated tasks, trying to minimize the amount of code you copy and
paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g.,
`run_exercise_3`).

Once you created a function that performs the given task, add it to the
`poisson-tester.cc` file, and make sure all the exercises are run through the
`gtest` executable, e.g., adding a test for each exercise, as in the following
snippet:

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have modified your Poisson code to
allow also non-homogeneous Neumann boundary conditions on different parts of
the domain, and you will have added some more options to the solver, enabling
usage of a direct solver, or of some more sofisticated preconditioners.

## Lab-06
### step-6

1.  See documentation of step-6 at
    <https://www.dealii.org/current/doxygen/deal.II/step_6.html>

2. Add the parameters
   
    - `Mapping degree`
    - `Marking strategy`
    - `Estimator type`
    - `Coarsening and refinement factors`
   
where `Mapping degree` controls the degree of the mapping used in the code,
`Marking strategy` is a choice between `global|fixed_fraction|fixed_number`,
`Estimator type` is a choice between `exact|kelly|residual`, and
`Coarsening and refinement factors` is a `std::pair<double, double>` containing
the arguments to pass to the `GridRefinement::refine_and_coarsen_fixed_*`
functions

1. Make sure all `FEValues` classes use a mapping with the correct order, and
make sure you use the correct mapping in the output as well (if you have a
recent Paraview)

4. Add a `Vector<float` field `error_per_cell` to `Poisson`, to be filled by
the method `estimate`

5. Add a method `residual_error_estimator` that computes the residual error estimator, using `FEInterfaceValues` and `FEValues`

6. Add a method `estimate` to the `Poisson` class, to compute the H1 seminorm
of the difference between the exact and computed solution if `Estimator type`
is `exact`, calls `KellyErrorEstimator<dim>::estimate` if `Estimator type` is
`kelly`, and calls `residual_error_estimator` if `Estimator type` is `residual`

7. Add to the convergence tables also the estimator you computed. This should be
identical to the `H1` semi-norm in the case where `Estimator type` is `exact`

8. Taking the exact solution computed on the L-shaped domain, compute the rate
at which the adaptive finite element method converges in terms of the number of
degrees of freedom using the three estimators above

9. Set zero boundary conditions, forcing term equal to four, exact solution
equal to `-x^2-y^2+1` and solve the problem on a circle with center in the
origin and radius one, for various finite element and mapping degrees. What do
you observe when the mapping degree does not match the finite element degree?
How do you explain this?

10. Create a test that reproduces exactly the behaviour of `step-6`, using only
your parameter file