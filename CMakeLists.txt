cmake_minimum_required(VERSION 3.24)

project(TyWRFCore VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

option(TYWRF_ENABLE_OPENMP "Enable OpenMP for CPU reference kernels when available" ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)

add_library(tywrf_core STATIC
  src/core/version.cpp
  src/dynamics/dynamics_loop.cpp
  src/io/netcdf_schema.cpp
)
add_library(TyWRF::Core ALIAS tywrf_core)

target_compile_features(tywrf_core PUBLIC cxx_std_20)
target_include_directories(tywrf_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_definitions(tywrf_core
  PUBLIC
    TYWRF_VERSION="${PROJECT_VERSION}"
)
target_link_libraries(tywrf_core PUBLIC PkgConfig::NETCDF)

if(TYWRF_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(tywrf_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(tywrf_core PUBLIC TYWRF_HAS_OPENMP=1)
  else()
    target_compile_definitions(tywrf_core PUBLIC TYWRF_HAS_OPENMP=0)
  endif()
else()
  target_compile_definitions(tywrf_core PUBLIC TYWRF_HAS_OPENMP=0)
endif()

set_target_properties(tywrf_core PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

if(BUILD_TESTING)
  add_executable(tywrf_smoke tests/unit/smoke_test.cpp)
  target_link_libraries(tywrf_smoke PRIVATE TyWRF::Core)
  add_test(NAME tywrf_smoke COMMAND tywrf_smoke)

  add_executable(tywrf_netcdf_schema_test tests/unit/netcdf_schema_test.cpp)
  target_link_libraries(tywrf_netcdf_schema_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_netcdf_schema_test COMMAND tywrf_netcdf_schema_test)

  add_executable(tywrf_dynamics_loop_test tests/unit/dynamics_loop_test.cpp)
  target_link_libraries(tywrf_dynamics_loop_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_dynamics_loop_test COMMAND tywrf_dynamics_loop_test)

  add_executable(tywrf_field_layout_test tests/unit/field_layout_test.cpp)
  target_link_libraries(tywrf_field_layout_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_field_layout_test COMMAND tywrf_field_layout_test)
endif()
