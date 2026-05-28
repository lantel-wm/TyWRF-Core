cmake_minimum_required(VERSION 3.24)

project(TyWRFCore VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

option(TYWRF_ENABLE_OPENMP "Enable OpenMP for CPU reference kernels when available" ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)

add_library(tywrf_core STATIC
  src/core/version.cpp
  src/dynamics/base_state.cpp
  src/dynamics/base_state_provider.cpp
  src/dynamics/cycle_schedule.cpp
  src/dynamics/dynamics_loop.cpp
  src/dynamics/pressure_refresh.cpp
  src/dynamics/pressure_refresh_hook.cpp
  src/dynamics/pressure_refresh_staging.cpp
  src/dynamics/state_stepper.cpp
  src/dynamics/tendency.cpp
  src/io/forcing_frames.cpp
  src/io/forcing_apply.cpp
  src/io/forcing_io.cpp
  src/io/netcdf_schema.cpp
  src/io/pressure_refresh_io.cpp
  src/io/wrf_state_io.cpp
  src/nest/base_state_exchange_adapter.cpp
  src/nest/base_state_exchange.cpp
  src/nest/base_state_source_staging.cpp
  src/nest/nest_interface.cpp
  src/nest/parent_child_interpolation.cpp
  src/nest/static_fields.cpp
  src/nest/state_exchange.cpp
  src/nest/state_remap.cpp
  src/nest/time_level_sync.cpp
  src/physics_bridge/staging.cpp
)
add_library(TyWRF::Core ALIAS tywrf_core)

target_compile_features(tywrf_core PUBLIC cxx_std_20)
target_include_directories(tywrf_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bindings/wrf_physics>
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

add_executable(tywrf_skeleton_cycle src/tools/skeleton_cycle.cpp)
target_link_libraries(tywrf_skeleton_cycle PRIVATE TyWRF::Core)
set_target_properties(tywrf_skeleton_cycle PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

add_executable(tywrf_tendency_cycle src/tools/tendency_cycle.cpp)
target_link_libraries(tywrf_tendency_cycle PRIVATE TyWRF::Core)
set_target_properties(tywrf_tendency_cycle PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

add_executable(tywrf_selected_field_cycle src/tools/selected_field_cycle.cpp)
target_link_libraries(tywrf_selected_field_cycle PRIVATE TyWRF::Core)
set_target_properties(tywrf_selected_field_cycle PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

add_executable(tywrf_pressure_refresh_probe src/tools/pressure_refresh_probe.cpp)
target_link_libraries(tywrf_pressure_refresh_probe PRIVATE TyWRF::Core)
set_target_properties(tywrf_pressure_refresh_probe PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

add_executable(tywrf_base_state_probe src/tools/base_state_probe.cpp)
target_link_libraries(tywrf_base_state_probe PRIVATE TyWRF::Core)
set_target_properties(tywrf_base_state_probe PROPERTIES
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

  add_executable(tywrf_cycle_schedule_test tests/unit/cycle_schedule_test.cpp)
  target_link_libraries(tywrf_cycle_schedule_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_cycle_schedule_test COMMAND tywrf_cycle_schedule_test)

  add_executable(tywrf_tendency_test tests/unit/tendency_test.cpp)
  target_link_libraries(tywrf_tendency_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_tendency_test COMMAND tywrf_tendency_test)

  add_executable(tywrf_state_stepper_test tests/unit/state_stepper_test.cpp)
  target_link_libraries(tywrf_state_stepper_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_state_stepper_test COMMAND tywrf_state_stepper_test)

  add_executable(tywrf_field_layout_test tests/unit/field_layout_test.cpp)
  target_link_libraries(tywrf_field_layout_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_field_layout_test COMMAND tywrf_field_layout_test)

  add_executable(tywrf_nest_interface_test tests/unit/nest_interface_test.cpp)
  target_link_libraries(tywrf_nest_interface_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_nest_interface_test COMMAND tywrf_nest_interface_test)

  add_executable(tywrf_state_remap_test tests/unit/state_remap_test.cpp)
  target_link_libraries(tywrf_state_remap_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_state_remap_test COMMAND tywrf_state_remap_test)

  add_executable(tywrf_state_exchange_test tests/unit/state_exchange_test.cpp)
  target_link_libraries(tywrf_state_exchange_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_state_exchange_test COMMAND tywrf_state_exchange_test)

  add_executable(tywrf_base_state_exchange_test
    tests/unit/base_state_exchange_test.cpp)
  target_link_libraries(tywrf_base_state_exchange_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_base_state_exchange_test
    COMMAND tywrf_base_state_exchange_test)

  add_executable(tywrf_base_state_exchange_adapter_test
    tests/unit/base_state_exchange_adapter_test.cpp)
  target_link_libraries(
    tywrf_base_state_exchange_adapter_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_base_state_exchange_adapter_test
    COMMAND tywrf_base_state_exchange_adapter_test)

  add_executable(tywrf_base_state_source_staging_test
    tests/unit/base_state_source_staging_test.cpp)
  target_link_libraries(
    tywrf_base_state_source_staging_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_base_state_source_staging_test
    COMMAND tywrf_base_state_source_staging_test)

  add_executable(tywrf_parent_child_interpolation_test
    tests/unit/parent_child_interpolation_test.cpp)
  target_link_libraries(tywrf_parent_child_interpolation_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_parent_child_interpolation_test
    COMMAND tywrf_parent_child_interpolation_test)

  add_executable(tywrf_static_fields_test tests/unit/static_fields_test.cpp)
  target_link_libraries(tywrf_static_fields_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_static_fields_test COMMAND tywrf_static_fields_test)

  add_executable(tywrf_time_level_sync_test tests/unit/time_level_sync_test.cpp)
  target_link_libraries(tywrf_time_level_sync_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_time_level_sync_test COMMAND tywrf_time_level_sync_test)

  add_executable(tywrf_pressure_refresh_test tests/unit/pressure_refresh_test.cpp)
  target_link_libraries(tywrf_pressure_refresh_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_pressure_refresh_test COMMAND tywrf_pressure_refresh_test)

  add_executable(tywrf_pressure_refresh_staging_test
    tests/unit/pressure_refresh_staging_test.cpp)
  target_link_libraries(tywrf_pressure_refresh_staging_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_pressure_refresh_staging_test
    COMMAND tywrf_pressure_refresh_staging_test)

  add_executable(tywrf_pressure_refresh_hook_test
    tests/unit/pressure_refresh_hook_test.cpp)
  target_link_libraries(tywrf_pressure_refresh_hook_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_pressure_refresh_hook_test
    COMMAND tywrf_pressure_refresh_hook_test)

  add_executable(tywrf_base_state_test tests/unit/base_state_test.cpp)
  target_link_libraries(tywrf_base_state_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_base_state_test COMMAND tywrf_base_state_test)

  add_executable(tywrf_base_state_provider_test
    tests/unit/base_state_provider_test.cpp)
  target_link_libraries(tywrf_base_state_provider_test PRIVATE TyWRF::Core)
  add_test(
    NAME tywrf_base_state_provider_test
    COMMAND tywrf_base_state_provider_test)

  add_executable(tywrf_pressure_refresh_io_test tests/unit/pressure_refresh_io_test.cpp)
  target_link_libraries(tywrf_pressure_refresh_io_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_pressure_refresh_io_test COMMAND tywrf_pressure_refresh_io_test)

  add_executable(tywrf_wrf_state_io_test tests/unit/wrf_state_io_test.cpp)
  target_link_libraries(tywrf_wrf_state_io_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_wrf_state_io_test COMMAND tywrf_wrf_state_io_test)

  add_executable(tywrf_forcing_io_test tests/unit/forcing_io_test.cpp)
  target_link_libraries(tywrf_forcing_io_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_forcing_io_test COMMAND tywrf_forcing_io_test)

  add_executable(tywrf_forcing_frames_test tests/unit/forcing_frames_test.cpp)
  target_link_libraries(tywrf_forcing_frames_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_forcing_frames_test COMMAND tywrf_forcing_frames_test)

  add_executable(tywrf_forcing_apply_test tests/unit/forcing_apply_test.cpp)
  target_link_libraries(tywrf_forcing_apply_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_forcing_apply_test COMMAND tywrf_forcing_apply_test)

  add_executable(tywrf_physics_bridge_test tests/unit/physics_bridge_test.cpp)
  target_link_libraries(tywrf_physics_bridge_test PRIVATE TyWRF::Core)
  add_test(NAME tywrf_physics_bridge_test COMMAND tywrf_physics_bridge_test)

  add_executable(tywrf_skeleton_cycle_tool_test tests/unit/skeleton_cycle_tool_test.cpp)
  target_link_libraries(tywrf_skeleton_cycle_tool_test PRIVATE PkgConfig::NETCDF)
  add_test(
    NAME tywrf_skeleton_cycle_tool_test
    COMMAND tywrf_skeleton_cycle_tool_test $<TARGET_FILE:tywrf_skeleton_cycle>
  )

  add_executable(tywrf_tendency_cycle_tool_test tests/unit/tendency_cycle_tool_test.cpp)
  target_link_libraries(tywrf_tendency_cycle_tool_test PRIVATE PkgConfig::NETCDF)
  add_test(
    NAME tywrf_tendency_cycle_tool_test
    COMMAND tywrf_tendency_cycle_tool_test $<TARGET_FILE:tywrf_tendency_cycle>
  )

  add_executable(tywrf_selected_field_cycle_tool_test
    tests/unit/selected_field_cycle_tool_test.cpp)
  target_link_libraries(tywrf_selected_field_cycle_tool_test PRIVATE PkgConfig::NETCDF)
  add_test(
    NAME tywrf_selected_field_cycle_tool_test
    COMMAND tywrf_selected_field_cycle_tool_test $<TARGET_FILE:tywrf_selected_field_cycle>
  )

  add_executable(tywrf_pressure_refresh_probe_tool_test
    tests/unit/pressure_refresh_probe_tool_test.cpp)
  target_link_libraries(tywrf_pressure_refresh_probe_tool_test PRIVATE PkgConfig::NETCDF)
  add_test(
    NAME tywrf_pressure_refresh_probe_tool_test
    COMMAND tywrf_pressure_refresh_probe_tool_test $<TARGET_FILE:tywrf_pressure_refresh_probe>
  )

  add_executable(tywrf_base_state_probe_tool_test
    tests/unit/base_state_probe_tool_test.cpp)
  target_link_libraries(tywrf_base_state_probe_tool_test PRIVATE PkgConfig::NETCDF)
  add_test(
    NAME tywrf_base_state_probe_tool_test
    COMMAND tywrf_base_state_probe_tool_test $<TARGET_FILE:tywrf_base_state_probe>
  )
endif()
