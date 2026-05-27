#include "tywrf/physics_bridge/staging.hpp"

#include <cstdint>

namespace tywrf::physics_bridge {
namespace {

constexpr std::int32_t kRequiredFieldCount = 25;
constexpr std::int32_t kRequiredSidecarV2FieldCount = 23;

struct ValidationResult {
  Status status = Status::stub_validated;
  std::int32_t failing_field = TYWRF_PHYSICS_FIELD_NONE;
  std::int32_t validated_field_count = 0;
};

void set_diagnostics(
    TywrfPhysicsDiagnostics* diagnostics,
    const ValidationResult result) noexcept {
  if (diagnostics == nullptr) {
    return;
  }

  diagnostics->status = static_cast<std::int32_t>(result.status);
  diagnostics->failing_field = result.failing_field;
  diagnostics->validated_field_count = result.validated_field_count;
  diagnostics->executed_physics = 0;
}

[[nodiscard]] constexpr bool is_supported_domain(const std::int32_t domain) noexcept {
  return domain == TYWRF_PHYSICS_DOMAIN_D01 || domain == TYWRF_PHYSICS_DOMAIN_D02;
}

[[nodiscard]] constexpr double abs_diff(const double lhs, const double rhs) noexcept {
  return lhs > rhs ? lhs - rhs : rhs - lhs;
}

[[nodiscard]] constexpr bool near_equal(const double lhs, const double rhs) noexcept {
  return abs_diff(lhs, rhs) <= 1.0e-9;
}

[[nodiscard]] ValidationResult failure(
    const Status status,
    const std::int32_t failing_field,
    const std::int32_t validated_field_count) noexcept {
  return {status, failing_field, validated_field_count};
}

[[nodiscard]] ValidationResult wrapper_unavailable(
    const std::int32_t validated_field_count) noexcept {
  return {Status::wrapper_unavailable, TYWRF_PHYSICS_FIELD_NONE, validated_field_count};
}

[[nodiscard]] ValidationResult validate_grid(const TywrfPhysicsGridMetadata& grid) noexcept {
  if (!is_supported_domain(grid.domain_id)) {
    return failure(Status::unsupported_domain, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (grid.mass_nx <= 0 || grid.mass_ny <= 0 || grid.mass_nz <= 0 ||
      grid.full_nz != grid.mass_nz + 1) {
    return failure(Status::invalid_dimensions, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  const auto expected_dx =
      static_cast<double>(grid.domain_id == TYWRF_PHYSICS_DOMAIN_D01 ? 10'000 : 2'000);
  const auto expected_dt =
      static_cast<double>(grid.domain_id == TYWRF_PHYSICS_DOMAIN_D01 ? 40 : 8);

  if (!near_equal(grid.dx_m, expected_dx) || !near_equal(grid.dy_m, expected_dx) ||
      !near_equal(grid.dt_s, expected_dt) || grid.step_index < 0 ||
      !near_equal(grid.end_seconds - grid.start_seconds, grid.dt_s)) {
    return failure(Status::invalid_timing, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  return {};
}

[[nodiscard]] ValidationResult validate_suite(
    const TywrfPhysicsSuiteConfig& suite,
    const std::int32_t domain) noexcept {
  const auto expected_cu = domain == TYWRF_PHYSICS_DOMAIN_D01 ? 1 : 0;

  if (suite.mp_physics != 8 || suite.cu_physics != expected_cu ||
      suite.ra_lw_physics != 4 || suite.ra_sw_physics != 4 ||
      suite.bl_pbl_physics != 1 || suite.sf_sfclay_physics != 1 ||
      suite.sf_surface_physics != 2 || suite.sf_ocean_physics != 1 ||
      suite.isftcflx != 2 || suite.num_moist_species != 8 ||
      suite.num_soil_layers != 4) {
    return failure(Status::unsupported_suite, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  return {};
}

[[nodiscard]] constexpr std::int32_t active_nx(
    const TywrfPhysicsField2D& field) noexcept {
  return field.nx - field.halo_i_lower - field.halo_i_upper;
}

[[nodiscard]] constexpr std::int32_t active_ny(
    const TywrfPhysicsField2D& field) noexcept {
  return field.ny - field.halo_j_lower - field.halo_j_upper;
}

[[nodiscard]] constexpr std::int32_t active_nx(
    const TywrfPhysicsField3D& field) noexcept {
  return field.nx - field.halo_i_lower - field.halo_i_upper;
}

[[nodiscard]] constexpr std::int32_t active_ny(
    const TywrfPhysicsField3D& field) noexcept {
  return field.ny - field.halo_j_lower - field.halo_j_upper;
}

[[nodiscard]] constexpr std::int32_t active_nz(
    const TywrfPhysicsField3D& field) noexcept {
  return field.nz - field.halo_k_lower - field.halo_k_upper;
}

[[nodiscard]] constexpr bool has_valid_element_size(const std::int32_t element_bytes) noexcept {
  return element_bytes == 4 || element_bytes == 8;
}

[[nodiscard]] ValidationResult validate_element_size(
    const std::int32_t element_bytes,
    const std::int32_t field_id,
    std::int32_t& expected_element_bytes,
    const std::int32_t validated_field_count) noexcept {
  if (!has_valid_element_size(element_bytes)) {
    return failure(Status::invalid_element_size, field_id, validated_field_count);
  }

  if (expected_element_bytes == 0) {
    expected_element_bytes = element_bytes;
    return {};
  }

  if (expected_element_bytes != element_bytes) {
    return failure(Status::invalid_element_size, field_id, validated_field_count);
  }

  return {};
}

[[nodiscard]] ValidationResult validate_field_2d(
    const TywrfPhysicsField2D& field,
    const std::int32_t field_id,
    const std::int32_t expected_active_nx,
    const std::int32_t expected_active_ny,
    std::int32_t& expected_element_bytes,
    std::int32_t& validated_field_count) noexcept {
  if (field.data == nullptr) {
    return failure(Status::missing_required_field, field_id, validated_field_count);
  }

  const auto element_result = validate_element_size(
      field.element_bytes, field_id, expected_element_bytes, validated_field_count);
  if (element_result.status != Status::stub_validated) {
    return element_result;
  }

  if (field.nx <= 0 || field.ny <= 0 || field.halo_i_lower < 0 ||
      field.halo_i_upper < 0 || field.halo_j_lower < 0 || field.halo_j_upper < 0 ||
      active_nx(field) != expected_active_nx || active_ny(field) != expected_active_ny) {
    return failure(Status::invalid_dimensions, field_id, validated_field_count);
  }

  if (field.stride_i != 1 || field.stride_j != field.nx) {
    return failure(Status::invalid_strides, field_id, validated_field_count);
  }

  ++validated_field_count;
  return {};
}

[[nodiscard]] ValidationResult validate_field_3d(
    const TywrfPhysicsField3D& field,
    const std::int32_t field_id,
    const std::int32_t expected_active_nx,
    const std::int32_t expected_active_ny,
    const std::int32_t expected_active_nz,
    std::int32_t& expected_element_bytes,
    std::int32_t& validated_field_count) noexcept {
  if (field.data == nullptr) {
    return failure(Status::missing_required_field, field_id, validated_field_count);
  }

  const auto element_result = validate_element_size(
      field.element_bytes, field_id, expected_element_bytes, validated_field_count);
  if (element_result.status != Status::stub_validated) {
    return element_result;
  }

  if (field.nx <= 0 || field.ny <= 0 || field.nz <= 0 || field.halo_i_lower < 0 ||
      field.halo_i_upper < 0 || field.halo_j_lower < 0 || field.halo_j_upper < 0 ||
      field.halo_k_lower < 0 || field.halo_k_upper < 0 ||
      active_nx(field) != expected_active_nx || active_ny(field) != expected_active_ny ||
      active_nz(field) != expected_active_nz) {
    return failure(Status::invalid_dimensions, field_id, validated_field_count);
  }

  if (field.stride_i != 1 || field.stride_k != field.nx ||
      field.stride_j != field.nx * field.nz) {
    return failure(Status::invalid_strides, field_id, validated_field_count);
  }

  ++validated_field_count;
  return {};
}

[[nodiscard]] ValidationResult validate_impl(const TywrfPhysicsStaging* staging) noexcept {
  if (staging == nullptr) {
    return failure(Status::null_argument, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (staging->abi_version != TYWRF_PHYSICS_ABI_VERSION) {
    return failure(Status::unsupported_abi, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  auto result = validate_grid(staging->grid);
  if (result.status != Status::stub_validated) {
    return result;
  }

  result = validate_suite(staging->suite, staging->grid.domain_id);
  if (result.status != Status::stub_validated) {
    return result;
  }

  const auto mass_nx = staging->grid.mass_nx;
  const auto mass_ny = staging->grid.mass_ny;
  const auto mass_nz = staging->grid.mass_nz;
  const auto full_nz = staging->grid.full_nz;

  std::int32_t expected_element_bytes = 0;
  std::int32_t validated_field_count = 0;

  result = validate_field_3d(
      staging->u, TYWRF_PHYSICS_FIELD_U, mass_nx + 1, mass_ny, mass_nz,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  result = validate_field_3d(
      staging->v, TYWRF_PHYSICS_FIELD_V, mass_nx, mass_ny + 1, mass_nz,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  result = validate_field_3d(
      staging->w, TYWRF_PHYSICS_FIELD_W, mass_nx, mass_ny, full_nz,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  result = validate_field_3d(
      staging->ph, TYWRF_PHYSICS_FIELD_PH, mass_nx, mass_ny, full_nz,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  result = validate_field_3d(
      staging->phb, TYWRF_PHYSICS_FIELD_PHB, mass_nx, mass_ny, full_nz,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }

  const TywrfPhysicsField3D* mass_fields[] = {
      &staging->t,      &staging->p,      &staging->pb,     &staging->qvapor,
      &staging->qcloud, &staging->qrain,  &staging->qice,   &staging->qsnow,
      &staging->qgraup, &staging->qnice,  &staging->qnrain,
  };
  const std::int32_t mass_field_ids[] = {
      TYWRF_PHYSICS_FIELD_T,      TYWRF_PHYSICS_FIELD_P,
      TYWRF_PHYSICS_FIELD_PB,     TYWRF_PHYSICS_FIELD_QVAPOR,
      TYWRF_PHYSICS_FIELD_QCLOUD, TYWRF_PHYSICS_FIELD_QRAIN,
      TYWRF_PHYSICS_FIELD_QICE,   TYWRF_PHYSICS_FIELD_QSNOW,
      TYWRF_PHYSICS_FIELD_QGRAUP, TYWRF_PHYSICS_FIELD_QNICE,
      TYWRF_PHYSICS_FIELD_QNRAIN,
  };

  for (std::int32_t index = 0; index < 11; ++index) {
    result = validate_field_3d(
        *mass_fields[index], mass_field_ids[index], mass_nx, mass_ny, mass_nz,
        expected_element_bytes, validated_field_count);
    if (result.status != Status::stub_validated) {
      return result;
    }
  }

  const TywrfPhysicsField2D* surface_fields[] = {
      &staging->mu,   &staging->mub,   &staging->psfc,  &staging->u10,   &staging->v10,
      &staging->t2,   &staging->q2,    &staging->rainc, &staging->rainnc,
  };
  const std::int32_t surface_field_ids[] = {
      TYWRF_PHYSICS_FIELD_MU,    TYWRF_PHYSICS_FIELD_MUB,
      TYWRF_PHYSICS_FIELD_PSFC,  TYWRF_PHYSICS_FIELD_U10,
      TYWRF_PHYSICS_FIELD_V10,   TYWRF_PHYSICS_FIELD_T2,
      TYWRF_PHYSICS_FIELD_Q2,    TYWRF_PHYSICS_FIELD_RAINC,
      TYWRF_PHYSICS_FIELD_RAINNC,
  };

  for (std::int32_t index = 0; index < 9; ++index) {
    result = validate_field_2d(
        *surface_fields[index], surface_field_ids[index], mass_nx, mass_ny,
        expected_element_bytes, validated_field_count);
    if (result.status != Status::stub_validated) {
      return result;
    }
  }

  return {Status::stub_validated, TYWRF_PHYSICS_FIELD_NONE, kRequiredFieldCount};
}

[[nodiscard]] ValidationResult validate_optional_field_2d(
    const TywrfPhysicsField2D& field,
    const std::int32_t field_id,
    const std::int32_t expected_active_nx,
    const std::int32_t expected_active_ny,
    std::int32_t& expected_element_bytes,
    std::int32_t& validated_field_count) noexcept {
  if (field.data == nullptr) {
    return {};
  }

  return validate_field_2d(
      field, field_id, expected_active_nx, expected_active_ny, expected_element_bytes,
      validated_field_count);
}

struct SidecarBlocksV2 {
  const TywrfPhysicsDriverContextV2* driver = nullptr;
  const TywrfPhysicsDerivedStateV2* derived = nullptr;
  const TywrfPhysicsStaticMaskV2* static_mask = nullptr;
  const TywrfPhysicsSfclaySurfaceV2* sfclay = nullptr;
};

[[nodiscard]] bool has_capability(
    const TywrfPhysicsBlockHeader& header,
    const std::uint64_t capability) noexcept {
  return (header.capabilities & capability) != 0U;
}

[[nodiscard]] ValidationResult validate_v2_header(
    const TywrfPhysicsBlockHeader* header,
    const std::uint32_t min_struct_size) noexcept {
  if (header == nullptr) {
    return failure(Status::null_argument, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (!tywrf_physics_block_header_has_min_size(
          header, TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE)) {
    return failure(Status::invalid_dimensions, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (!tywrf_physics_block_header_has_abi(header, TYWRF_PHYSICS_ABI_VERSION_V2)) {
    return failure(Status::unsupported_abi, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (!tywrf_physics_block_header_has_min_size(header, min_struct_size)) {
    return failure(Status::invalid_dimensions, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  return {};
}

[[nodiscard]] ValidationResult collect_sidecar_blocks_v2(
    const TywrfPhysicsBlockHeader* sidecars,
    SidecarBlocksV2& blocks) noexcept {
  if (sidecars == nullptr) {
    return failure(Status::missing_required_field, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  const auto* header = sidecars;
  while (header != nullptr) {
    auto result = validate_v2_header(header, TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE);
    if (result.status != Status::stub_validated) {
      return result;
    }

    if (has_capability(*header, TYWRF_PHYSICS_CAPABILITY_DRIVER_CONTEXT)) {
      result = validate_v2_header(header, TYWRF_PHYSICS_DRIVER_CONTEXT_V2_STRUCT_SIZE);
      if (result.status != Status::stub_validated) {
        return result;
      }
      blocks.driver = reinterpret_cast<const TywrfPhysicsDriverContextV2*>(header);
    } else if (has_capability(*header, TYWRF_PHYSICS_CAPABILITY_DERIVED_STATE)) {
      result = validate_v2_header(header, TYWRF_PHYSICS_DERIVED_STATE_V2_STRUCT_SIZE);
      if (result.status != Status::stub_validated) {
        return result;
      }
      blocks.derived = reinterpret_cast<const TywrfPhysicsDerivedStateV2*>(header);
    } else if (has_capability(*header, TYWRF_PHYSICS_CAPABILITY_STATIC_MASK)) {
      result = validate_v2_header(header, TYWRF_PHYSICS_STATIC_MASK_V2_STRUCT_SIZE);
      if (result.status != Status::stub_validated) {
        return result;
      }
      blocks.static_mask = reinterpret_cast<const TywrfPhysicsStaticMaskV2*>(header);
    } else if (has_capability(*header, TYWRF_PHYSICS_CAPABILITY_SURFACE_STATE)) {
      result = validate_v2_header(header, TYWRF_PHYSICS_SFCLAY_SURFACE_V2_STRUCT_SIZE);
      if (result.status != Status::stub_validated) {
        return result;
      }
      blocks.sfclay = reinterpret_cast<const TywrfPhysicsSfclaySurfaceV2*>(header);
    }

    header = static_cast<const TywrfPhysicsBlockHeader*>(header->next);
  }

  if (blocks.driver == nullptr || blocks.derived == nullptr ||
      blocks.static_mask == nullptr || blocks.sfclay == nullptr) {
    return failure(Status::missing_required_field, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  return {};
}

[[nodiscard]] ValidationResult validate_driver_context_v2(
    const TywrfPhysicsStaging& staging,
    const TywrfPhysicsDriverContextV2& driver) noexcept {
  const TywrfPhysicsGridMetadata grid{
      driver.domain_id,
      driver.mass_nx,
      driver.mass_ny,
      driver.mass_nz,
      driver.full_nz,
      driver.dx_m,
      driver.dy_m,
      driver.dt_s,
      driver.step_index,
      driver.start_seconds,
      driver.end_seconds,
  };

  auto result = validate_grid(grid);
  if (result.status != Status::stub_validated) {
    return result;
  }

  result = validate_suite(driver.suite, driver.domain_id);
  if (result.status != Status::stub_validated) {
    return result;
  }

  if (driver.domain_id != staging.grid.domain_id || driver.mass_nx != staging.grid.mass_nx ||
      driver.mass_ny != staging.grid.mass_ny || driver.mass_nz != staging.grid.mass_nz ||
      driver.full_nz != staging.grid.full_nz || !near_equal(driver.dx_m, staging.grid.dx_m) ||
      !near_equal(driver.dy_m, staging.grid.dy_m) || !near_equal(driver.dt_s, staging.grid.dt_s)) {
    return failure(Status::invalid_dimensions, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (driver.ids > driver.ide || driver.jds > driver.jde || driver.kds > driver.kde ||
      driver.ims > driver.ime || driver.jms > driver.jme || driver.kms > driver.kme ||
      driver.its > driver.ite || driver.jts > driver.jte || driver.kts > driver.kte ||
      driver.its < driver.ims || driver.ite > driver.ime || driver.jts < driver.jms ||
      driver.jte > driver.jme || driver.kts < driver.kms || driver.kte > driver.kme) {
    return failure(Status::invalid_dimensions, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  if (driver.enable_sfclay == 0 || driver.constants == nullptr) {
    return failure(Status::missing_required_field, TYWRF_PHYSICS_FIELD_NONE, 0);
  }

  return {};
}

[[nodiscard]] ValidationResult validate_derived_state_v2(
    const TywrfPhysicsDerivedStateV2& derived,
    const TywrfPhysicsDriverContextV2& driver,
    std::int32_t& validated_field_count) noexcept {
  const TywrfPhysicsField3D* fields[] = {
      &derived.u_phy,
      &derived.v_phy,
      &derived.t_phy,
      &derived.qv_curr,
      &derived.p_phy,
      &derived.dz8w,
  };
  const std::int32_t field_ids[] = {
      TYWRF_PHYSICS_FIELD_U_PHY,
      TYWRF_PHYSICS_FIELD_V_PHY,
      TYWRF_PHYSICS_FIELD_T_PHY,
      TYWRF_PHYSICS_FIELD_QV_CURR,
      TYWRF_PHYSICS_FIELD_P_PHY,
      TYWRF_PHYSICS_FIELD_DZ8W,
  };

  std::int32_t expected_element_bytes = 0;
  for (std::int32_t index = 0; index < 6; ++index) {
    const auto result = validate_field_3d(
        *fields[index], field_ids[index], driver.mass_nx, driver.mass_ny,
        driver.mass_nz, expected_element_bytes, validated_field_count);
    if (result.status != Status::stub_validated) {
      return result;
    }
  }

  return {};
}

[[nodiscard]] ValidationResult validate_static_mask_v2(
    const TywrfPhysicsStaticMaskV2& static_mask,
    const TywrfPhysicsDriverContextV2& driver,
    std::int32_t& validated_field_count) noexcept {
  std::int32_t expected_element_bytes = 0;

  auto result = validate_field_2d(
      static_mask.xland, TYWRF_PHYSICS_FIELD_XLAND, driver.mass_nx, driver.mass_ny,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  expected_element_bytes = 0;
  result = validate_field_2d(
      static_mask.lakemask, TYWRF_PHYSICS_FIELD_LAKEMASK, driver.mass_nx,
      driver.mass_ny, expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  expected_element_bytes = 0;
  result = validate_field_2d(
      static_mask.lu_index, TYWRF_PHYSICS_FIELD_LU_INDEX, driver.mass_nx,
      driver.mass_ny, expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  expected_element_bytes = 0;
  result = validate_field_2d(
      static_mask.hgt, TYWRF_PHYSICS_FIELD_HGT, driver.mass_nx, driver.mass_ny,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }

  expected_element_bytes = 0;
  result = validate_optional_field_2d(
      static_mask.xlat, TYWRF_PHYSICS_FIELD_XLAT, driver.mass_nx, driver.mass_ny,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }

  expected_element_bytes = 0;
  return validate_optional_field_2d(
      static_mask.xlong, TYWRF_PHYSICS_FIELD_XLONG, driver.mass_nx, driver.mass_ny,
      expected_element_bytes, validated_field_count);
}

[[nodiscard]] ValidationResult validate_sfclay_surface_v2(
    const TywrfPhysicsSfclaySurfaceV2& sfclay,
    const TywrfPhysicsDriverContextV2& driver,
    std::int32_t& validated_field_count) noexcept {
  const TywrfPhysicsField2D* fields[] = {
      &sfclay.psfc, &sfclay.tsk, &sfclay.znt, &sfclay.ust, &sfclay.pblh,
      &sfclay.hfx,  &sfclay.qfx, &sfclay.lh,  &sfclay.u10, &sfclay.v10,
      &sfclay.t2,   &sfclay.q2,  &sfclay.th2,
  };
  const std::int32_t field_ids[] = {
      TYWRF_PHYSICS_FIELD_PSFC, TYWRF_PHYSICS_FIELD_TSK,
      TYWRF_PHYSICS_FIELD_ZNT,  TYWRF_PHYSICS_FIELD_UST,
      TYWRF_PHYSICS_FIELD_PBLH, TYWRF_PHYSICS_FIELD_HFX,
      TYWRF_PHYSICS_FIELD_QFX,  TYWRF_PHYSICS_FIELD_LH,
      TYWRF_PHYSICS_FIELD_U10,  TYWRF_PHYSICS_FIELD_V10,
      TYWRF_PHYSICS_FIELD_T2,   TYWRF_PHYSICS_FIELD_Q2,
      TYWRF_PHYSICS_FIELD_TH2,
  };

  std::int32_t expected_element_bytes = 0;
  for (std::int32_t index = 0; index < 13; ++index) {
    const auto result = validate_field_2d(
        *fields[index], field_ids[index], driver.mass_nx, driver.mass_ny,
        expected_element_bytes, validated_field_count);
    if (result.status != Status::stub_validated) {
      return result;
    }
  }

  auto result = validate_optional_field_2d(
      sfclay.scratch_1, TYWRF_PHYSICS_FIELD_NONE, driver.mass_nx, driver.mass_ny,
      expected_element_bytes, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }

  return validate_optional_field_2d(
      sfclay.scratch_2, TYWRF_PHYSICS_FIELD_NONE, driver.mass_nx, driver.mass_ny,
      expected_element_bytes, validated_field_count);
}

[[nodiscard]] ValidationResult validate_sidecar_impl_v2(
    const TywrfPhysicsStaging* staging,
    const TywrfPhysicsBlockHeader* sidecars) noexcept {
  auto result = validate_impl(staging);
  if (result.status != Status::stub_validated) {
    return result;
  }

  SidecarBlocksV2 blocks{};
  result = collect_sidecar_blocks_v2(sidecars, blocks);
  if (result.status != Status::stub_validated) {
    return result;
  }

  result = validate_driver_context_v2(*staging, *blocks.driver);
  if (result.status != Status::stub_validated) {
    return result;
  }

  std::int32_t validated_field_count = 0;
  result = validate_derived_state_v2(*blocks.derived, *blocks.driver, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  result = validate_static_mask_v2(*blocks.static_mask, *blocks.driver, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }
  result = validate_sfclay_surface_v2(*blocks.sfclay, *blocks.driver, validated_field_count);
  if (result.status != Status::stub_validated) {
    return result;
  }

  return {
      Status::stub_validated, TYWRF_PHYSICS_FIELD_NONE,
      validated_field_count >= kRequiredSidecarV2FieldCount ? validated_field_count
                                                            : kRequiredSidecarV2FieldCount};
}

}  // namespace

Status validate_staging(
    const TywrfPhysicsStaging& staging,
    TywrfPhysicsDiagnostics* diagnostics) noexcept {
  const auto result = validate_impl(&staging);
  set_diagnostics(diagnostics, result);
  return result.status;
}

Status validate_sidecar_v2(
    const TywrfPhysicsStaging& staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics) noexcept {
  const auto result = validate_sidecar_impl_v2(&staging, sidecars);
  set_diagnostics(diagnostics, result);
  return result.status;
}

Status run_stub_bridge(
    const TywrfPhysicsStaging& staging,
    TywrfPhysicsDiagnostics* diagnostics) noexcept {
  return static_cast<Status>(tywrf_wrf_physics_step(&staging, diagnostics));
}

Status run_stub_bridge_ex(
    const TywrfPhysicsStaging& staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics) noexcept {
  return static_cast<Status>(tywrf_wrf_physics_step_ex(&staging, sidecars, diagnostics));
}

std::string_view status_name(const Status status) noexcept {
  return tywrf_physics_status_name(static_cast<std::int32_t>(status));
}

}  // namespace tywrf::physics_bridge

extern "C" int32_t tywrf_wrf_physics_step(
    const TywrfPhysicsStaging* staging,
    TywrfPhysicsDiagnostics* diagnostics) {
  const auto result = tywrf::physics_bridge::validate_impl(staging);
  tywrf::physics_bridge::set_diagnostics(diagnostics, result);
  return static_cast<std::int32_t>(result.status);
}

extern "C" int32_t tywrf_physics_validate_sidecar_v2(
    const TywrfPhysicsStaging* staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics) {
  const auto result = tywrf::physics_bridge::validate_sidecar_impl_v2(staging, sidecars);
  tywrf::physics_bridge::set_diagnostics(diagnostics, result);
  return static_cast<std::int32_t>(result.status);
}

extern "C" int32_t tywrf_wrf_physics_step_ex(
    const TywrfPhysicsStaging* staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics) {
  const auto result = tywrf::physics_bridge::validate_sidecar_impl_v2(staging, sidecars);
  if (result.status != tywrf::physics_bridge::Status::stub_validated) {
    tywrf::physics_bridge::set_diagnostics(diagnostics, result);
    return static_cast<std::int32_t>(result.status);
  }

  const auto unavailable =
      tywrf::physics_bridge::wrapper_unavailable(result.validated_field_count);
  tywrf::physics_bridge::set_diagnostics(diagnostics, unavailable);
  return static_cast<std::int32_t>(unavailable.status);
}

extern "C" const char* tywrf_physics_status_name(const int32_t status) {
  switch (status) {
    case TYWRF_PHYSICS_STATUS_STUB_VALIDATED:
      return "stub_validated";
    case TYWRF_PHYSICS_STATUS_NULL_ARGUMENT:
      return "null_argument";
    case TYWRF_PHYSICS_STATUS_UNSUPPORTED_ABI:
      return "unsupported_abi";
    case TYWRF_PHYSICS_STATUS_UNSUPPORTED_DOMAIN:
      return "unsupported_domain";
    case TYWRF_PHYSICS_STATUS_UNSUPPORTED_SUITE:
      return "unsupported_suite";
    case TYWRF_PHYSICS_STATUS_INVALID_TIMING:
      return "invalid_timing";
    case TYWRF_PHYSICS_STATUS_INVALID_DIMENSIONS:
      return "invalid_dimensions";
    case TYWRF_PHYSICS_STATUS_INVALID_STRIDES:
      return "invalid_strides";
    case TYWRF_PHYSICS_STATUS_INVALID_ELEMENT_SIZE:
      return "invalid_element_size";
    case TYWRF_PHYSICS_STATUS_MISSING_REQUIRED_FIELD:
      return "missing_required_field";
    case TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE:
      return "wrapper_unavailable";
    default:
      return "unknown";
  }
}
