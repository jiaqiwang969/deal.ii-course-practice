#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <gtest/gtest.h>
// 只对 proc 0 进行监听
int
main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  testing::InitGoogleTest(&argc, argv);

  ::testing::TestEventListeners &listeners =
    ::testing::UnitTest::GetInstance()->listeners();

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) != 0)
    {
      delete listeners.Release(listeners.default_result_printer());
    }
  return RUN_ALL_TESTS();
}