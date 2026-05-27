#include "tywrf/io/pressure_refresh_io.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct CliOptions {
  std::filesystem::path input;
  std::int32_t mass_nx = 0;
  std::int32_t mass_ny = 0;
  std::int32_t mass_nz = 0;
  std::int32_t full_nz = 0;
  std::int32_t halo = 0;
  std::size_t time_index = 0;
  bool pretty = false;
};

[[nodiscard]] std::string usage() {
  return "usage: tywrf_pressure_refresh_probe --input PATH --mass-nx N --mass-ny N "
         "--mass-nz N --full-nz N [--halo N] [--time-index N] [--pretty]";
}

[[nodiscard]] std::int32_t parse_i32(const std::string_view value, const std::string_view name) {
  std::size_t parsed = 0;
  int parsed_value = 0;
  try {
    parsed_value = std::stoi(std::string(value), &parsed);
  } catch (const std::exception&) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  if (parsed != value.size()) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  return static_cast<std::int32_t>(parsed_value);
}

[[nodiscard]] std::size_t parse_size(const std::string_view value, const std::string_view name) {
  std::size_t parsed = 0;
  unsigned long long parsed_value = 0;
  try {
    parsed_value = std::stoull(std::string(value), &parsed);
  } catch (const std::exception&) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  if (parsed != value.size()) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  return static_cast<std::size_t>(parsed_value);
}

[[nodiscard]] CliOptions parse_args(const int argc, char** argv) {
  CliOptions options;
  for (int index = 1; index < argc; ++index) {
    const std::string_view arg(argv[index]);
    const auto require_value = [&](const std::string_view name) -> std::string_view {
      if (index + 1 >= argc) {
        throw std::runtime_error("missing value for " + std::string(name));
      }
      ++index;
      return argv[index];
    };

    if (arg == "--input") {
      options.input = std::filesystem::path(require_value(arg));
    } else if (arg == "--mass-nx") {
      options.mass_nx = parse_i32(require_value(arg), arg);
    } else if (arg == "--mass-ny") {
      options.mass_ny = parse_i32(require_value(arg), arg);
    } else if (arg == "--mass-nz") {
      options.mass_nz = parse_i32(require_value(arg), arg);
    } else if (arg == "--full-nz") {
      options.full_nz = parse_i32(require_value(arg), arg);
    } else if (arg == "--halo") {
      options.halo = parse_i32(require_value(arg), arg);
    } else if (arg == "--time-index") {
      options.time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--pretty") {
      options.pretty = true;
    } else if (arg == "--help" || arg == "-h") {
      std::cout << usage() << '\n';
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + std::string(arg));
    }
  }

  if (options.input.empty()) {
    throw std::runtime_error("missing required --input");
  }
  if (options.mass_nx <= 0 || options.mass_ny <= 0 || options.mass_nz <= 0 ||
      options.full_nz <= 0) {
    throw std::runtime_error("--mass-nx, --mass-ny, --mass-nz, and --full-nz must be positive");
  }
  if (options.halo < 0) {
    throw std::runtime_error("--halo must be non-negative");
  }
  return options;
}

[[nodiscard]] const char* p_top_source_name(
    const tywrf::io::PressureRefreshPTopSource source) noexcept {
  switch (source) {
    case tywrf::io::PressureRefreshPTopSource::missing:
      return "missing";
    case tywrf::io::PressureRefreshPTopSource::global_attribute:
      return "global_attribute";
    case tywrf::io::PressureRefreshPTopSource::scalar_variable:
      return "scalar_variable";
    case tywrf::io::PressureRefreshPTopSource::time_variable:
      return "time_variable";
  }
  return "unknown";
}

[[nodiscard]] std::string json_escape(const std::string_view value) {
  std::ostringstream escaped;
  for (const char c : value) {
    switch (c) {
      case '\\':
        escaped << "\\\\";
        break;
      case '"':
        escaped << "\\\"";
        break;
      case '\n':
        escaped << "\\n";
        break;
      case '\r':
        escaped << "\\r";
        break;
      case '\t':
        escaped << "\\t";
        break;
      default:
        escaped << c;
        break;
    }
  }
  return escaped.str();
}

void write_indent(std::ostream& output, const bool pretty, const int level) {
  if (!pretty) {
    return;
  }
  output << '\n';
  for (int index = 0; index < level; ++index) {
    output << "  ";
  }
}

void write_string_array(
    std::ostream& output,
    const std::vector<std::string>& values,
    const bool pretty,
    const int level) {
  output << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    write_indent(output, pretty, level + 1);
    output << '"' << json_escape(values[index]) << '"';
  }
  if (!values.empty()) {
    write_indent(output, pretty, level);
  }
  output << ']';
}

void write_shape3(
    std::ostream& output,
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t nz) {
  output << '[' << nx << ',' << ny << ',' << nz << ']';
}

void write_key(std::ostream& output, const std::string_view key, const bool pretty) {
  output << '"' << key << "\":";
  if (pretty) {
    output << ' ';
  }
}

template <typename Writer>
void write_field(
    std::ostream& output,
    const std::string_view key,
    const bool pretty,
    bool& first,
    Writer writer) {
  if (!first) {
    output << ',';
  }
  first = false;
  write_indent(output, pretty, 1);
  write_key(output, key, pretty);
  writer();
}

[[nodiscard]] std::string status_for_report(
    const tywrf::io::KrosaPressureRefreshInputReport& report) {
  return report.ok() ? "ok" : "missing";
}

void write_report_json(
    std::ostream& output,
    const tywrf::io::KrosaPressureRefreshInputReport& report,
    const bool pretty) {
  bool first = true;
  output << '{';
  write_field(output, "status", pretty, first, [&] {
    output << '"' << status_for_report(report) << '"';
  });
  write_field(output, "input_file", pretty, first, [&] {
    output << '"' << json_escape(report.source_path.string()) << '"';
  });
  write_field(output, "time_index", pretty, first, [&] { output << report.time_index; });
  write_field(output, "expected_counts", pretty, first, [&] {
    output << '{';
    bool counts_first = true;
    write_field(output, "mass_levels", pretty, counts_first, [&] {
      output << report.expected_mass_level_count;
    });
    write_field(output, "full_levels", pretty, counts_first, [&] {
      output << report.expected_full_level_count;
    });
    write_field(output, "alb_points", pretty, counts_first, [&] {
      output << report.expected_alb_point_count;
    });
    write_indent(output, pretty, 1);
    output << '}';
  });
  write_field(output, "actual_counts", pretty, first, [&] {
    output << '{';
    bool counts_first = true;
    write_field(output, "c3f", pretty, counts_first, [&] { output << report.c3f_count; });
    write_field(output, "c4f", pretty, counts_first, [&] { output << report.c4f_count; });
    write_field(output, "c3h", pretty, counts_first, [&] { output << report.c3h_count; });
    write_field(output, "c4h", pretty, counts_first, [&] { output << report.c4h_count; });
    write_field(output, "alb_points", pretty, counts_first, [&] {
      output << report.alb_point_count;
    });
    write_indent(output, pretty, 1);
    output << '}';
  });
  write_field(output, "expected_shapes", pretty, first, [&] {
    output << '{';
    bool shapes_first = true;
    write_field(output, "alb", pretty, shapes_first, [&] {
      write_shape3(
          output,
          report.expected_alb_nx,
          report.expected_alb_ny,
          report.expected_alb_nz);
    });
    write_field(output, "mass_levels", pretty, shapes_first, [&] {
      output << '[' << report.expected_mass_level_count << ']';
    });
    write_field(output, "full_levels", pretty, shapes_first, [&] {
      output << '[' << report.expected_full_level_count << ']';
    });
    write_indent(output, pretty, 1);
    output << '}';
  });
  write_field(output, "actual_shapes", pretty, first, [&] {
    output << '{';
    bool shapes_first = true;
    write_field(output, "alb", pretty, shapes_first, [&] {
      write_shape3(output, report.alb_nx, report.alb_ny, report.alb_nz);
    });
    write_field(output, "c3f", pretty, shapes_first, [&] {
      output << '[' << report.c3f_count << ']';
    });
    write_field(output, "c4f", pretty, shapes_first, [&] {
      output << '[' << report.c4f_count << ']';
    });
    write_field(output, "c3h", pretty, shapes_first, [&] {
      output << '[' << report.c3h_count << ']';
    });
    write_field(output, "c4h", pretty, shapes_first, [&] {
      output << '[' << report.c4h_count << ']';
    });
    write_indent(output, pretty, 1);
    output << '}';
  });
  write_field(output, "missing_names", pretty, first, [&] {
    write_string_array(output, report.missing_names, pretty, 1);
  });
  write_field(output, "p_top_present", pretty, first, [&] {
    output << (report.p_top_present() ? "true" : "false");
  });
  write_field(output, "p_top_source", pretty, first, [&] {
    output << '"' << p_top_source_name(report.p_top_source) << '"';
  });
  write_field(output, "alb_available", pretty, first, [&] {
    output << (report.alb_point_count == report.expected_alb_point_count ? "true" : "false");
  });
  write_field(output, "alb_loaded", pretty, first, [&] {
    output << (report.alb_loaded ? "true" : "false");
  });
  write_field(output, "diagnostic_only", pretty, first, [&] { output << "true"; });
  write_field(output, "calls_pressure_refresh_compute", pretty, first, [&] {
    output << "false";
  });
  write_field(output, "reads_or_writes_p", pretty, first, [&] { output << "false"; });
  write_indent(output, pretty, 0);
  output << '}';
  if (pretty) {
    output << '\n';
  }
}

void write_error_json(
    std::ostream& output,
    const CliOptions& options,
    const std::string_view message,
    const bool pretty) {
  bool first = true;
  output << '{';
  write_field(output, "status", pretty, first, [&] { output << "\"error\""; });
  write_field(output, "input_file", pretty, first, [&] {
    output << '"' << json_escape(options.input.string()) << '"';
  });
  write_field(output, "time_index", pretty, first, [&] { output << options.time_index; });
  write_field(output, "error", pretty, first, [&] {
    output << '"' << json_escape(message) << '"';
  });
  write_field(output, "diagnostic_only", pretty, first, [&] { output << "true"; });
  write_field(output, "calls_pressure_refresh_compute", pretty, first, [&] {
    output << "false";
  });
  write_field(output, "reads_or_writes_p", pretty, first, [&] { output << "false"; });
  write_indent(output, pretty, 0);
  output << '}';
  if (pretty) {
    output << '\n';
  }
}

[[nodiscard]] tywrf::Grid make_grid(const CliOptions& options) {
  return tywrf::Grid(
      {options.mass_nx,
       options.mass_ny,
       options.mass_nz,
       options.full_nz,
       tywrf::uniform_halo_3d(options.halo)});
}

}  // namespace

int main(const int argc, char** argv) {
  CliOptions options;
  try {
    options = parse_args(argc, argv);
    const auto grid = make_grid(options);
    tywrf::FieldStorage3D<float> alb(grid.mass_layout());
    const auto result = tywrf::io::read_krosa_pressure_refresh_inputs(
        options.input,
        grid,
        alb,
        {.time_index = options.time_index});
    write_report_json(std::cout, result.report, options.pretty);
    return result.ok() ? 0 : 1;
  } catch (const std::exception& error) {
    if (options.input.empty()) {
      std::cerr << usage() << '\n';
      std::cerr << "pressure-refresh probe error: " << error.what() << '\n';
      return 2;
    }
    write_error_json(std::cout, options, error.what(), options.pretty);
    std::cerr << "pressure-refresh probe error: " << error.what() << '\n';
    return 1;
  }
}
