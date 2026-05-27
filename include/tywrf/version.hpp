#pragma once

#include <string_view>

namespace tywrf {

[[nodiscard]] std::string_view project_name() noexcept;
[[nodiscard]] std::string_view project_version() noexcept;
[[nodiscard]] bool openmp_enabled() noexcept;

}  // namespace tywrf
