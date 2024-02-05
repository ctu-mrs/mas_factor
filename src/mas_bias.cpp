#include "mas_factor/mas_bias.h"

#include <gtsam/geometry/Point3.h>
#include <iostream>

namespace gtsam {

namespace mas_bias {

std::ostream& operator<<(std::ostream& os, const ConstantBias& bias) {
  os << "lin_acc = " << bias.linAcc().transpose();
  os << "ang_acc = " << bias.angAcc().transpose();
  return os;
}

/// print with optional string
void ConstantBias::print(const std::string& s) const {
  std::cout << s << *this << std::endl;
}

} // namespace mas_bias

} // namespace gtsam

