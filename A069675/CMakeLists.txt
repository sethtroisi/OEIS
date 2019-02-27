cmake_minimum_required(VERSION 2.8.12)
project(A069675)

# g++ -O3 -g A069675_sieve.cpp -lgmp -lgmpxx -lprofiler --std=c++11
set(CMAKE_CXX_FLAGS "-O2 -std=c++11 ${CMAKE_CXX_FLAGS}")

add_subdirectory(abseil-cpp)

include(FindOpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-g ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(sieve A069675_sieve.cpp)
target_link_libraries(sieve gmp)
target_link_libraries(sieve gmpxx)
target_link_libraries(sieve absl::container absl::hash)
#target_link_libraries(sieve profiler)

#######################################

add_executable(tester A069675_tester.cpp)
target_link_libraries(tester gmp)
target_link_libraries(tester gmpxx)
#target_link_libraries(tester profiler)
