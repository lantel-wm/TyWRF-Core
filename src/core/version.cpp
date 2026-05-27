#include "tywrf/version.hpp"

#ifndef TYWRF_VERSION
#define TYWRF_VERSION "0.1.0"
#endif

#ifndef TYWRF_HAS_OPENMP
#define TYWRF_HAS_OPENMP 0
#endif

namespace tywrf {

std::string_view project_name() noexcept {
  return "TyWRF-Core";
}

std::string_view project_version() noexcept {
  return TYWRF_VERSION;
}

bool openmp_enabled() noexcept {
  return TYWRF_HAS_OPENMP != 0;
}

}  // namespace tywrf
