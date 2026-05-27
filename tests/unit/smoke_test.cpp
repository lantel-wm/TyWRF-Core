#include "tywrf/version.hpp"

#include <cassert>
#include <iostream>
#include <string_view>

int main() {
  assert(tywrf::project_name() == std::string_view{"TyWRF-Core"});
  assert(!tywrf::project_version().empty());

  std::cout << tywrf::project_name() << " " << tywrf::project_version()
            << " openmp=" << (tywrf::openmp_enabled() ? "on" : "off")
            << '\n';
  return 0;
}
