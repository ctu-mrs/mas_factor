#pragma once

#include <gtsam/base/OptionalJacobian.h>
#include <gtsam/base/VectorSpace.h>
#include <iosfwd>
#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
#include <boost/serialization/nvp.hpp>
#endif

namespace gtsam {

namespace mas_bias {

class GTSAM_EXPORT ConstantBias {
private:
  Vector3 bias_lin_acc_; ///< The units for stddev are σ = m/s² or m √Hz/s²
  Vector3 bias_ang_acc_; ///< The units for stddev are σ = rad/s² or rad √Hz/s²

public:
  /// dimension of the variable - used to autodetect sizes
  static const size_t dimension = 6;

  /// @name Standard Constructors
  /// @{

  ConstantBias() :
      bias_lin_acc_(0.0, 0.0, 0.0), bias_ang_acc_(0.0, 0.0, 0.0) {
  }

  ConstantBias(const Vector3& bias_lin_acc, const Vector3& bias_ang_acc) :
      bias_lin_acc_(bias_lin_acc), bias_ang_acc_(bias_ang_acc) {
  }

  explicit ConstantBias(const Vector6& v) :
      bias_lin_acc_(v.head<3>()), bias_ang_acc_(v.tail<3>()) {
  }

  /// @}

  /** return the linear and angular acceleration biases in a single vector */
  Vector6 vector() const {
    Vector6 v;
    v << bias_lin_acc_, bias_ang_acc_;
    return v;
  }

  /** get linear acceleration bias */
  const Vector3& linAcc() const {
    return bias_lin_acc_;
  }

  /** get angular acceleration bias */
  const Vector3& angAcc() const {
    return bias_ang_acc_;
  }

  /** Correct a linear acceleration measurement using this bias model, and optionally compute Jacobians */
  Vector3 correctLinAcc(const Vector3& measurement,
                               OptionalJacobian<3, 6> H1 = {},
                               OptionalJacobian<3, 3> H2 = {}) const {
    if (H1) (*H1) << -I_3x3, Z_3x3;
    if (H2) (*H2) << I_3x3;
    return measurement - bias_lin_acc_;
  }

  /** Correct an angular acceleration measurement using this bias model, and optionally compute Jacobians */
  Vector3 correctAngAcc(const Vector3& measurement,
                           OptionalJacobian<3, 6> H1 = {},
                           OptionalJacobian<3, 3> H2 = {}) const {
    if (H1) (*H1) << Z_3x3, -I_3x3;
    if (H2) (*H2) << I_3x3;
    return measurement - bias_ang_acc_;
  }

  /// @name Testable
  /// @{

  /// ostream operator
  GTSAM_EXPORT friend std::ostream& operator<<(std::ostream& os,
                                               const ConstantBias& bias);

  /// print with optional string
  void print(const std::string& s = "") const;

  /** equality up to tolerance */
  inline bool equals(const ConstantBias& expected, double tol = 1e-5) const {
    return equal_with_abs_tol(bias_lin_acc_, expected.bias_lin_acc_, tol)
        && equal_with_abs_tol(bias_ang_acc_, expected.bias_ang_acc_, tol);
  }

  /// @}
  /// @name Group
  /// @{

  /** identity for group operation */
  static ConstantBias Identity() {
    return ConstantBias();
  }

  /** inverse */
  inline ConstantBias operator-() const {
    return ConstantBias(-bias_lin_acc_, -bias_ang_acc_);
  }

  /** addition of vector on right */
  ConstantBias operator+(const Vector6& v) const {
    return ConstantBias(bias_lin_acc_ + v.head<3>(), bias_ang_acc_ + v.tail<3>());
  }

  /** addition */
  ConstantBias operator+(const ConstantBias& b) const {
    return ConstantBias(bias_lin_acc_ + b.bias_lin_acc_, bias_ang_acc_ + b.bias_ang_acc_);
  }

  /** subtraction */
  ConstantBias operator-(const ConstantBias& b) const {
    return ConstantBias(bias_lin_acc_ - b.bias_lin_acc_, bias_ang_acc_ - b.bias_ang_acc_);
  }

  /// @}

private:

  /// @name Advanced Interface
  /// @{

#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar & BOOST_SERIALIZATION_NVP(bias_lin_acc_);
    ar & BOOST_SERIALIZATION_NVP(bias_ang_acc_);
  }
#endif


public:
  GTSAM_MAKE_ALIGNED_OPERATOR_NEW
  /// @}

}; // ConstantBias class
} // namespace mas_bias

template<>
struct traits<mas_bias::ConstantBias> : public internal::VectorSpace<
    mas_bias::ConstantBias> {
};

} // namespace gtsam

