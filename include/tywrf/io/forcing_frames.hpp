#pragma once

#include "tywrf/io/forcing_io.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace tywrf::io {

inline constexpr std::int64_t kKrosaForcingIntervalSeconds = 21'600;
inline constexpr std::size_t kMaxForcingFieldRank = 4;

class ForcingFrameError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

struct ForcingFrame {
  std::size_t record_index = 0;
  std::int64_t start_seconds = 0;
  std::int64_t end_seconds = 0;
  std::string start_time;
  std::string end_time;
};

struct ForcingFrameWeights {
  double old_weight = 1.0;
  double new_weight = 0.0;
};

struct ForcingFrameSelectorConfig {
  std::size_t time_count = 0;
  std::int64_t interval_seconds = kKrosaForcingIntervalSeconds;
  std::vector<std::string> time_strings;
};

class ForcingFrameSelector {
 public:
  explicit ForcingFrameSelector(ForcingFrameSelectorConfig config);

  [[nodiscard]] std::size_t time_count() const noexcept;
  [[nodiscard]] std::int64_t interval_seconds() const noexcept;
  [[nodiscard]] std::span<const std::string> time_strings() const noexcept;

  [[nodiscard]] ForcingFrame frame_for_cycle_index(std::size_t cycle_index) const;
  [[nodiscard]] ForcingFrame frame_for_start_seconds(std::int64_t start_seconds) const;
  [[nodiscard]] ForcingFrame frame_containing_model_seconds(
      std::int64_t model_seconds) const;
  [[nodiscard]] ForcingFrameWeights weights_for(
      const ForcingFrame& frame,
      std::int64_t model_seconds) const;

 private:
  ForcingFrameSelectorConfig config_;
};

struct ForcingFieldLayout {
  std::size_t rank = 0;
  std::array<std::size_t, kMaxForcingFieldRank> extents{};
  std::array<std::size_t, kMaxForcingFieldRank> strides{};
  std::size_t value_count = 0;
};

struct ForcingFieldView {
  float* data = nullptr;
  ForcingFieldLayout layout{};

  [[nodiscard]] std::size_t index(std::span<const std::size_t> indices) const;
  [[nodiscard]] float& at(std::span<const std::size_t> indices) const;
};

struct ForcingFieldBuffer {
  ForcingFieldLayout layout{};
  std::vector<float> values;

  [[nodiscard]] ForcingFieldView view() noexcept;
};

struct ForcingStagedSlice {
  std::string variable_name;
  std::size_t record_index = 0;
  ForcingFieldBuffer field;
};

struct CanonicalForcingLayout {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::size_t value_count = 0;

  [[nodiscard]] std::size_t index(
      std::size_t i,
      std::size_t j,
      std::size_t k) const;
};

struct PackedForcingField {
  CanonicalForcingLayout layout{};
  std::vector<float> values;

  [[nodiscard]] float at(std::size_t i, std::size_t j, std::size_t k) const;
};

struct SpectralNudgingVariableMapping {
  std::string canonical_name;
  std::string old_variable_name;
  std::string new_variable_name;
};

[[nodiscard]] ForcingFrameSelector make_krosa_forcing_frame_selector(
    const KrosaForcingReader& reader,
    std::int64_t interval_seconds = kKrosaForcingIntervalSeconds);

[[nodiscard]] std::string add_seconds_to_wrf_timestamp(
    std::string_view timestamp,
    std::int64_t offset_seconds);

[[nodiscard]] ForcingFieldLayout make_forcing_field_layout(
    std::span<const std::size_t> shape);

[[nodiscard]] ForcingFieldBuffer make_forcing_field_buffer(
    const ForcingFieldLayout& layout);

[[nodiscard]] ForcingStagedSlice stage_forcing_time_slice(
    const ForcingTimeSlice& slice);

[[nodiscard]] CanonicalForcingLayout make_canonical_forcing_layout(
    std::size_t nx,
    std::size_t ny,
    std::size_t nz);

[[nodiscard]] PackedForcingField pack_fdda_3d_raw_to_canonical(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape);

[[nodiscard]] PackedForcingField pack_boundary_x_side_raw_to_canonical(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape);

[[nodiscard]] PackedForcingField pack_boundary_y_side_raw_to_canonical(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape);

[[nodiscard]] std::vector<SpectralNudgingVariableMapping>
krosa_spectral_nudging_variable_mappings();

}  // namespace tywrf::io
