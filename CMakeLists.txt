cmake_minimum_required(VERSION 2.8)
project(chaospp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(/usr/local/include ${PROJECT_SOURCE_DIR}/dependencies/include ${PROJECT_SOURCE_DIR}/chaospp)
link_directories(/usr/local/lib)

##### Sample

add_executable(escape_time_uniform_sample examples/escape_time_uniform_sample.cpp)
target_link_libraries(escape_time_uniform_sample gmp mpfr)

add_executable(escape_time_wang_landau examples/escape_time_wang_landau.cpp)
target_link_libraries(escape_time_wang_landau gmp mpfr)

add_executable(print_phase_space examples/print_phase_space.cpp)
target_link_libraries(print_phase_space gmp mpfr)

add_executable(lyapunov_sm6_uniform examples/lyapunov_sm6_uniform.cpp)
target_link_libraries(lyapunov_sm6_uniform gmp mpfr)

##### Search

add_executable(stagger search/powerlaw.cpp)
target_link_libraries(stagger gmp mpfr)

add_executable(adaptive search/adaptive.cpp)
target_link_libraries(adaptive gmp mpfr)

add_executable(anisotropic search/anisotropic.cpp)
target_link_libraries(anisotropic gmp mpfr)

#### Sample

add_executable(lyapunov_tent sample/lyapunov_tent.cpp)
target_link_libraries(lyapunov_tent gmp mpfr)

add_executable(lyapunov_sm6 sample/lyapunov_sm6.cpp)
target_link_libraries(lyapunov_sm6 gmp mpfr)

add_executable(wl_sy sample/wl_sy.cpp)
target_link_libraries(wl_sy gmp mpfr)

add_executable(divergence_wl_sy sample/divergence_wl_sy.cpp)
target_link_libraries(divergence_wl_sy gmp mpfr)

add_executable(efficiency_chm sample/efficiency_chm.cpp)
target_link_libraries(efficiency_chm gmp mpfr)

#### Test Assumptions

add_executable(tstar_control test_assumptions/tstar_control.cpp)
target_link_libraries(tstar_control gmp mpfr)

add_executable(lambda_correlation test_assumptions/lambda_correlation.cpp)
target_link_libraries(lambda_correlation gmp mpfr)

add_executable(tstar_lambda_tent test_assumptions/tstar_lambda_tent.cpp)
target_link_libraries(tstar_lambda_tent gmp mpfr)

add_executable(tstar_te test_assumptions/tstar_te.cpp)
target_link_libraries(tstar_te gmp mpfr)

#### Test canonical

add_executable(test_te_ce test_canonical/escape_time_canonical_ensemble.cpp)
target_link_libraries(test_te_ce gmp mpfr)

#### Tests

option(build_tests "Build the project's unit tests" ON)
if(build_tests)
    enable_testing()
    set(GTEST_ROOT ${PROJECT_SOURCE_DIR}/dependencies)

    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})

    add_executable(runUnitTests test/main.cpp)
    target_link_libraries(runUnitTests ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} gmp mpfr)
    add_test(runUnitTests runUnitTests)
endif()
