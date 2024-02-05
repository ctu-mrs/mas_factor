#include "mas_params.h"

using namespace std;

namespace gtsam {

//------------------------------------------------------------------------------
void MasParams::print(const string& s) const {
  PreintegrationParams::print(s);
  cout << "mass:\n[\n" << mass << "\n]"
       << endl;
  cout << "motorConstant:\n[\n" << motorConstant << "\n]"
       << endl;
  cout << "momentConstant:\n[\n" << momentConstant << "\n]"
       << endl;
  cout << "numRotors:\n[\n" << numRotors << "\n]"
       << endl;
  cout << "rotorDirs:\n[\n" << rotorDirs << "\n]"
       << endl;
}

//------------------------------------------------------------------------------
bool MasParams::equals(const PreintegrationParams& other,
                                  double tol) const {
  auto e = dynamic_cast<const MasParams*>(&other);
  return e != nullptr && PreintegrationParams::equals(other, tol) &&
         equal_with_abs_tol(mass, e->mass,
                            tol) &&
         equal_with_abs_tol(motorConstant, e->motorConstant,
                            tol) &&
         equal_with_abs_tol(momentConstant, e->momentConstant,
                            tol) &&
         equal_with_abs_tol(numRotors, e->numRotors,
                            tol);
}

}  // namespace gtsam
