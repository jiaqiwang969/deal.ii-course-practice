#include <deal.II/base/utilities.h>

#include "poisson.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize init(argc, argv);

  std::string par_name = "";
  if (argc > 1)
    par_name = argv[1];

  deallog.depth_console(2);
  Poisson<2> laplace_problem;
  laplace_problem.initialize(par_name);
  laplace_problem.run();
  return 0;
}
