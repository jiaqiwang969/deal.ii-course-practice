CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

PROJECT(gtest)

FIND_PACKAGE(deal.II 9.2 REQUIRED
HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Fix warning on mac
SET(CMAKE_MACOSX_RPATH 1)

# Executable for the exerciese
ADD_EXECUTABLE(poisson
    source/base_problem.cc 
    source/base_block_problem.cc 
    source/poisson.cc 
    source/base_problem.cc 
    source/linear_elasticity.cc
    source/stokes.cc
    source/main.cc)


DEAL_II_SETUP_TARGET(poisson)

ADD_EXECUTABLE(linear_elasticity
    source/base_problem.cc 
    source/base_block_problem.cc 
    source/poisson.cc 
    source/base_problem.cc 
    source/linear_elasticity.cc
    source/stokes.cc
    source/main.cc)


DEAL_II_SETUP_TARGET(linear_elasticity)

ADD_EXECUTABLE(stokes
    source/base_problem.cc 
    source/base_block_problem.cc 
    source/poisson.cc 
    source/base_problem.cc 
    source/linear_elasticity.cc
    source/stokes.cc
    source/main.cc)


DEAL_II_SETUP_TARGET(stokes)



# # Library of the executable
# ADD_LIBRARY(poisson-lib SHARED $<TARGET_OBJECTS:poisson>)
# DEAL_II_SETUP_TARGET(poisson-lib)

# # Tester executable
# FIND_PACKAGE(GTest)
# FILE(GLOB test_files tests/*cc)
# ADD_EXECUTABLE(gtest ${test_files})
# TARGET_LINK_LIBRARIES(gtest ${GTEST_LIBRARY} poisson-lib)   
# DEAL_II_SETUP_TARGET(gtest)

INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS} ./include/)

# ENABLE_TESTING()
# GTEST_DISCOVER_TESTS(gtest)