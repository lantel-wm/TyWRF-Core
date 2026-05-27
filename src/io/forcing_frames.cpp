#include "tywrf/io/forcing_frames.hpp"

#include <chrono>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace tywrf::io {
namespace {

void require(const bool condition, const char* message) {
  if (!condition) {
    throw ForcingFrameError(message);
  }
}

[[nodiscard]] int parse_fixed_int(
    const std::string_view value,
    const std::size_t offset,
    const std::size_t count) {
  int result = 0;
  for (std::size_t pos = 0; pos < count; ++pos) {
    const char c = value.at(offset + pos);
    if (c < '0' || c > '9') {
      throw ForcingFrameError("WRF timestamp contains a non-digit field");
    }
    result = result * 10 + static_cast<int>(c - '0');
  }
  return result;
}

[[nodiscard]] std::chrono::sys_seconds parse_wrf_timestamp(
    const std::string_view timestamp) {
  if (timestamp.size() != 19 || timestamp.at(4) != '-' ||
      timestamp.at(7) != '-' || timestamp.at(10) != '_' ||
      timestamp.at(13) != ':' || timestamp.at(16) != ':') {
    throw ForcingFrameError("WRF timestamp must use YYYY-MM-DD_HH:MM:SS");
  }

  const auto year = parse_fixed_int(timestamp, 0, 4);
  const auto month = parse_fixed_int(timestamp, 5, 2);
  const auto day = parse_fixed_int(timestamp, 8, 2);
  const auto hour = parse_fixed_int(timestamp, 11, 2);
  const auto minute = parse_fixed_int(timestamp, 14, 2);
  const auto second = parse_fixed_int(timestamp, 17, 2);

  const std::chrono::year_month_day ymd{
      std::chrono::year{year},
      std::chrono::month{static_cast<unsigned>(month)},
      std::chrono::day{static_cast<unsigned>(day)}};
  if (!ymd.ok() || hour > 23 || minute > 59 || second > 59) {
    throw ForcingFrameError("WRF timestamp has an out-of-range date or time");
  }

  return std::chrono::sys_days{ymd} + std::chrono::hours{hour} +
         std::chrono::minutes{minute} + std::chrono::seconds{second};
}

[[nodiscard]] std::string format_wrf_timestamp(const std::chrono::sys_seconds timestamp) {
  const auto day_floor = std::chrono::floor<std::chrono::days>(timestamp);
  const std::chrono::year_month_day ymd{day_floor};
  const std::chrono::hh_mm_ss tod{timestamp - day_floor};

  std::ostringstream out;
  out << std::setfill('0') << std::setw(4) << static_cast<int>(ymd.year())
      << '-' << std::setw(2) << static_cast<unsigned>(ymd.month())
      << '-' << std::setw(2) << static_cast<unsigned>(ymd.day())
      << '_' << std::setw(2) << tod.hours().count()
      << ':' << std::setw(2) << tod.minutes().count()
      << ':' << std::setw(2) << tod.seconds().count();
  return out.str();
}

void validate_config(const ForcingFrameSelectorConfig& config) {
  require(config.time_count > 0, "forcing selector requires at least one Time record");
  require(config.interval_seconds > 0, "forcing interval must be positive");
  require(!config.time_strings.empty(),
          "forcing selector requires explicit start timestamp strings");
  require(config.time_strings.size() == config.time_count,
          "forcing timestamp count must match Time dimension length");
  for (const auto& timestamp : config.time_strings) {
    (void)parse_wrf_timestamp(timestamp);
  }
}

[[nodiscard]] std::size_t checked_record_index(
    const std::int64_t start_seconds,
    const std::int64_t interval_seconds,
    const std::size_t time_count) {
  if (start_seconds < 0) {
    throw ForcingFrameError("forcing frame start_seconds must be non-negative");
  }
  if (start_seconds % interval_seconds != 0) {
    throw ForcingFrameError("forcing frame start_seconds must align to the interval");
  }
  const auto record_index = static_cast<std::size_t>(start_seconds / interval_seconds);
  if (record_index >= time_count) {
    throw ForcingFrameError("forcing frame start_seconds is outside available records");
  }
  return record_index;
}

[[nodiscard]] std::size_t checked_product(std::span<const std::size_t> shape) {
  std::size_t product = 1;
  for (const auto extent : shape) {
    if (extent != 0 &&
        product > std::numeric_limits<std::size_t>::max() / extent) {
      throw ForcingFrameError("forcing field layout size overflows size_t");
    }
    product *= extent;
  }
  return product;
}

[[nodiscard]] ForcingFieldLayout validate_raw_pack_input(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape) {
  if (raw_shape.size() != 3) {
    throw ForcingFrameError("raw forcing pack requires rank 3 input");
  }
  for (const auto extent : raw_shape) {
    if (extent == 0) {
      throw ForcingFrameError("raw forcing pack requires non-empty extents");
    }
  }

  auto raw_layout = make_forcing_field_layout(raw_shape);
  if (raw_values.size() != raw_layout.value_count) {
    throw ForcingFrameError("raw forcing value count does not match input shape");
  }
  return raw_layout;
}

[[nodiscard]] std::size_t raw_index_3d(
    const ForcingFieldLayout& layout,
    const std::size_t dim0,
    const std::size_t dim1,
    const std::size_t dim2) {
  return dim0 * layout.strides[0] + dim1 * layout.strides[1] +
         dim2 * layout.strides[2];
}

}  // namespace

ForcingFrameSelector::ForcingFrameSelector(ForcingFrameSelectorConfig config)
    : config_(std::move(config)) {
  validate_config(config_);
}

std::size_t ForcingFrameSelector::time_count() const noexcept {
  return config_.time_count;
}

std::int64_t ForcingFrameSelector::interval_seconds() const noexcept {
  return config_.interval_seconds;
}

std::span<const std::string> ForcingFrameSelector::time_strings() const noexcept {
  return config_.time_strings;
}

ForcingFrame ForcingFrameSelector::frame_for_cycle_index(
    const std::size_t cycle_index) const {
  if (cycle_index >= config_.time_count) {
    throw ForcingFrameError("forcing cycle index is outside available records");
  }

  const auto start_seconds =
      static_cast<std::int64_t>(cycle_index) * config_.interval_seconds;
  ForcingFrame frame;
  frame.record_index = cycle_index;
  frame.start_seconds = start_seconds;
  frame.end_seconds = start_seconds + config_.interval_seconds;
  frame.start_time = config_.time_strings.at(cycle_index);
  frame.end_time =
      cycle_index + 1 < config_.time_strings.size()
          ? config_.time_strings.at(cycle_index + 1)
          : add_seconds_to_wrf_timestamp(frame.start_time, config_.interval_seconds);
  return frame;
}

ForcingFrame ForcingFrameSelector::frame_for_start_seconds(
    const std::int64_t start_seconds) const {
  return frame_for_cycle_index(
      checked_record_index(start_seconds, config_.interval_seconds, config_.time_count));
}

ForcingFrame ForcingFrameSelector::frame_containing_model_seconds(
    const std::int64_t model_seconds) const {
  if (model_seconds < 0) {
    throw ForcingFrameError("model_seconds must be non-negative");
  }

  const auto final_end_seconds =
      static_cast<std::int64_t>(config_.time_count) * config_.interval_seconds;
  if (model_seconds > final_end_seconds) {
    throw ForcingFrameError("model_seconds is outside available forcing frames");
  }

  std::size_t cycle_index = 0;
  if (model_seconds > 0) {
    cycle_index = static_cast<std::size_t>(
        (model_seconds - 1) / config_.interval_seconds);
  }
  return frame_for_cycle_index(cycle_index);
}

ForcingFrameWeights ForcingFrameSelector::weights_for(
    const ForcingFrame& frame,
    const std::int64_t model_seconds) const {
  if (model_seconds < frame.start_seconds || model_seconds > frame.end_seconds) {
    throw ForcingFrameError("model_seconds is outside the selected forcing frame");
  }
  if (frame.end_seconds <= frame.start_seconds) {
    throw ForcingFrameError("forcing frame has a non-positive duration");
  }

  const auto elapsed = static_cast<double>(model_seconds - frame.start_seconds);
  const auto duration = static_cast<double>(frame.end_seconds - frame.start_seconds);
  const double new_weight = elapsed / duration;
  return {1.0 - new_weight, new_weight};
}

ForcingFrameSelector make_krosa_forcing_frame_selector(
    const KrosaForcingReader& reader,
    const std::int64_t interval_seconds) {
  return ForcingFrameSelector(
      ForcingFrameSelectorConfig{reader.time_count(), interval_seconds,
                                 reader.read_time_strings()});
}

std::string add_seconds_to_wrf_timestamp(
    const std::string_view timestamp,
    const std::int64_t offset_seconds) {
  return format_wrf_timestamp(
      parse_wrf_timestamp(timestamp) + std::chrono::seconds{offset_seconds});
}

ForcingFieldLayout make_forcing_field_layout(std::span<const std::size_t> shape) {
  if (shape.size() > kMaxForcingFieldRank) {
    throw ForcingFrameError("forcing field rank exceeds supported staging rank");
  }

  ForcingFieldLayout layout;
  layout.rank = shape.size();
  for (std::size_t dim = 0; dim < shape.size(); ++dim) {
    layout.extents.at(dim) = shape[dim];
  }
  layout.value_count = checked_product(shape);
  if (layout.rank == 0) {
    return layout;
  }

  layout.strides.at(layout.rank - 1) = 1;
  for (std::size_t dim = layout.rank - 1; dim > 0; --dim) {
    layout.strides.at(dim - 1) = layout.strides.at(dim) * layout.extents.at(dim);
  }
  return layout;
}

ForcingFieldBuffer make_forcing_field_buffer(const ForcingFieldLayout& layout) {
  ForcingFieldBuffer buffer;
  buffer.layout = layout;
  buffer.values.assign(layout.value_count, 0.0F);
  return buffer;
}

std::size_t ForcingFieldView::index(std::span<const std::size_t> indices) const {
  if (indices.size() != layout.rank) {
    throw ForcingFrameError("forcing field index rank mismatch");
  }
  std::size_t linear_index = 0;
  for (std::size_t dim = 0; dim < layout.rank; ++dim) {
    if (indices[dim] >= layout.extents.at(dim)) {
      throw ForcingFrameError("forcing field index is outside layout extent");
    }
    linear_index += indices[dim] * layout.strides.at(dim);
  }
  return linear_index;
}

float& ForcingFieldView::at(std::span<const std::size_t> indices) const {
  return data[index(indices)];
}

ForcingFieldView ForcingFieldBuffer::view() noexcept {
  return {values.data(), layout};
}

std::size_t CanonicalForcingLayout::index(
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const {
  if (i >= nx || j >= ny || k >= nz) {
    throw ForcingFrameError("canonical forcing index is outside layout extent");
  }
  return ((j * nz) + k) * nx + i;
}

float PackedForcingField::at(
    const std::size_t i,
    const std::size_t j,
    const std::size_t k) const {
  const auto linear_index = layout.index(i, j, k);
  if (linear_index >= values.size()) {
    throw ForcingFrameError("packed forcing value count does not match layout");
  }
  return values[linear_index];
}

ForcingStagedSlice stage_forcing_time_slice(const ForcingTimeSlice& slice) {
  ForcingStagedSlice staged;
  staged.variable_name = slice.metadata.name;
  staged.record_index = slice.time_index;
  staged.field.layout = make_forcing_field_layout(slice.metadata.slice_shape);
  staged.field.values = slice.values;
  if (staged.field.values.size() != staged.field.layout.value_count) {
    throw ForcingFrameError("forcing time-slice value count does not match layout");
  }
  return staged;
}

CanonicalForcingLayout make_canonical_forcing_layout(
    const std::size_t nx,
    const std::size_t ny,
    const std::size_t nz) {
  if (nx == 0 || ny == 0 || nz == 0) {
    throw ForcingFrameError("canonical forcing layout requires non-empty extents");
  }
  const std::array<std::size_t, 3> shape{ny, nz, nx};
  return CanonicalForcingLayout{nx, ny, nz, checked_product(shape)};
}

PackedForcingField pack_fdda_3d_raw_to_canonical(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape) {
  const auto raw_layout = validate_raw_pack_input(raw_values, raw_shape);
  const auto nz = raw_shape[0];
  const auto ny = raw_shape[1];
  const auto nx = raw_shape[2];

  PackedForcingField packed;
  packed.layout = make_canonical_forcing_layout(nx, ny, nz);
  packed.values.assign(packed.layout.value_count, 0.0F);

  for (std::size_t k = 0; k < nz; ++k) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t i = 0; i < nx; ++i) {
        packed.values[packed.layout.index(i, j, k)] =
            raw_values[raw_index_3d(raw_layout, k, j, i)];
      }
    }
  }
  return packed;
}

PackedForcingField pack_boundary_x_side_raw_to_canonical(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape) {
  const auto raw_layout = validate_raw_pack_input(raw_values, raw_shape);
  const auto bdy_width = raw_shape[0];
  const auto nz = raw_shape[1];
  const auto ny = raw_shape[2];

  PackedForcingField packed;
  packed.layout = make_canonical_forcing_layout(bdy_width, ny, nz);
  packed.values.assign(packed.layout.value_count, 0.0F);

  for (std::size_t bdy = 0; bdy < bdy_width; ++bdy) {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
        packed.values[packed.layout.index(bdy, j, k)] =
            raw_values[raw_index_3d(raw_layout, bdy, k, j)];
      }
    }
  }
  return packed;
}

PackedForcingField pack_boundary_y_side_raw_to_canonical(
    std::span<const float> raw_values,
    std::span<const std::size_t> raw_shape) {
  const auto raw_layout = validate_raw_pack_input(raw_values, raw_shape);
  const auto bdy_width = raw_shape[0];
  const auto nz = raw_shape[1];
  const auto nx = raw_shape[2];

  PackedForcingField packed;
  packed.layout = make_canonical_forcing_layout(nx, bdy_width, nz);
  packed.values.assign(packed.layout.value_count, 0.0F);

  for (std::size_t bdy = 0; bdy < bdy_width; ++bdy) {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t i = 0; i < nx; ++i) {
        packed.values[packed.layout.index(i, bdy, k)] =
            raw_values[raw_index_3d(raw_layout, bdy, k, i)];
      }
    }
  }
  return packed;
}

std::vector<SpectralNudgingVariableMapping> krosa_spectral_nudging_variable_mappings() {
  return {
      {"U", "U_NDG_OLD", "U_NDG_NEW"},
      {"V", "V_NDG_OLD", "V_NDG_NEW"},
      {"T", "T_NDG_OLD", "T_NDG_NEW"},
      {"QVAPOR", "Q_NDG_OLD", "Q_NDG_NEW"},
      {"PH", "PH_NDG_OLD", "PH_NDG_NEW"},
      {"MU", "MU_NDG_OLD", "MU_NDG_NEW"},
  };
}

}  // namespace tywrf::io
