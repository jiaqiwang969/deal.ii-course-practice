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