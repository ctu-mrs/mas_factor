#pragma once

#include <gtsam/geometry/Pose3.h>
#include <gtsam/base/Vector.h>
#include <gtsam/base/Manifold.h>

namespace gtsam
{

/// Velocity is currently typedef'd to Vector3
typedef Vector3 Velocity3;

typedef Eigen::Matrix<double, 10, 10> Matrix_10x10;
typedef Eigen::Matrix<double, 3, 10>  Matrix_3x10;
typedef Eigen::Matrix<double, 12, 12> Matrix_12x12;
typedef Eigen::Matrix<double, 3, 12>  Matrix_3x12;

class GTSAM_EXPORT FullState {
private:
  Rot3      R_;  ///< Rotation nRb, rotates points/velocities in body to points/velocities in nav
  Point3    t_;  ///< position n_t, in nav frame
  Velocity3 v_;  ///< linear velocity n_v in nav frame
  Velocity3 w_;  ///< angular velocity n_w in nav frame

public:
  enum
  {
    dimension = 12
  };

  typedef std::pair<Point3, Velocity3> PositionAndVelocity;

  /// @name Constructors
  /// @{

  /// Default constructor
  FullState() : t_(0, 0, 0), v_(Vector3::Zero()) {
  }
  /// Construct from attitude, position, velocity
  FullState(const Rot3& R, const Point3& t, const Velocity3& v, const Velocity3& w) : R_(R), t_(t), v_(v), w_(w) {
  }
  /// Construct from pose and velocity
  FullState(const Pose3& pose, const Velocity3& v, const Velocity3& w) : R_(pose.rotation()), t_(pose.translation()), v_(v), w_(w) {
  }
  /// Construct from SO(3) and R^6
  FullState(const Matrix3& R, const Vector6& tv, const Velocity3& w) : R_(R), t_(tv.head<3>()), v_(tv.tail<3>()), w_(w) {
  }
  /// Named constructor with derivatives
  static FullState Create(const Rot3& R, const Point3& t, const Velocity3& v, const Velocity3& w, OptionalJacobian<12, 3> H1, OptionalJacobian<12, 3> H2,
                          OptionalJacobian<12, 3> H3, OptionalJacobian<12, 3> H4);
  /// Named constructor with derivatives
  static FullState FromPoseVelocity(const Pose3& pose, const Vector3& lin_vel, const Vector3& ang_vel, OptionalJacobian<12, 6> H1, OptionalJacobian<12, 3> H2,
                                    OptionalJacobian<12, 3> H3);

  /// @}
  /// @name Component Access
  /// @{

  const Rot3&      attitude(OptionalJacobian<3, 12> H = boost::none) const;
  const Point3&    position(OptionalJacobian<3, 12> H = boost::none) const;
  const Velocity3& linVelocity(OptionalJacobian<3, 12> H = boost::none) const;
  const Velocity3& angVelocity(OptionalJacobian<3, 12> H = boost::none) const;

  const Pose3 pose() const {
    return Pose3(attitude(), position());
  }

  /// @}
  /// @name Derived quantities
  /// @{

  /// Return rotation matrix. Induces computation in quaternion mode
  Matrix3 R() const {
    return R_.matrix();
  }
  /// Return quaternion. Induces computation in matrix mode
  Quaternion quaternion() const {
    return R_.toQuaternion();
  }
  /// Return position as Vector3
  Vector3 t() const {
    return t_;
  }
  /// Return linear velocity as Vector3. Computation-free.
  const Vector3& v() const {
    return v_;
  }

  /// Return angular velocity as Vector3. Computation-free.
  const Vector3& w() const {
    return w_;
  }
  // Return velocity in body frame
  Velocity3 bodyLinVelocity(OptionalJacobian<3, 12> H = boost::none) const;

  // Return velocity in body frame
  Velocity3 bodyAngVelocity(OptionalJacobian<3, 12> H = boost::none) const;

  /// Return matrix group representation, in MATLAB notation:
  /// nTb = [nRb 0 n_t; 0 nRb n_v; 0 0 1]
  /// With this embedding in GL(3), matrix product agrees with compose
  Matrix_10x10 matrix() const;

  /// @}
  /// @name Testable
  /// @{

  /// Output stream operator
  GTSAM_EXPORT
  friend std::ostream& operator<<(std::ostream& os, const FullState& state);

  /// print
  void print(const std::string& s = "") const;

  /// equals
  bool equals(const FullState& other, double tol = 1e-8) const;

  /// @}
  /// @name Manifold
  /// @{

  // Tangent space sugar.
  static Eigen::Block<Vector12, 3, 1> dR(Vector12& v) {
    return v.segment<3>(0);
  }
  static Eigen::Block<Vector12, 3, 1> dP(Vector12& v) {
    return v.segment<3>(3);
  }
  static Eigen::Block<Vector12, 3, 1> dV(Vector12& v) {
    return v.segment<3>(6);
  }
  static Eigen::Block<Vector12, 3, 1> dW(Vector12& v) {
    return v.segment<3>(9);
  }
  static Eigen::Block<const Vector12, 3, 1> dR(const Vector12& v) {
    return v.segment<3>(0);
  }
  static Eigen::Block<const Vector12, 3, 1> dP(const Vector12& v) {
    return v.segment<3>(3);
  }
  static Eigen::Block<const Vector12, 3, 1> dV(const Vector12& v) {
    return v.segment<3>(6);
  }
  static Eigen::Block<const Vector12, 3, 1> dW(const Vector12& v) {
    return v.segment<3>(9);
  }

  /// retract with optional derivatives
  FullState retract(const Vector12&          v,  //
                    OptionalJacobian<12, 12> H1 = boost::none, OptionalJacobian<12, 12> H2 = boost::none) const;

  /// localCoordinates with optional derivatives
  Vector12 localCoordinates(const FullState&         g,  //
                            OptionalJacobian<12, 12> H1 = boost::none, OptionalJacobian<12, 12> H2 = boost::none) const;

  /// @}
  /// @name Dynamics
  /// @{

  /// Integrate forward in time given angular acceleration and linear acceleration in body frame
  /// Uses second order integration for position, returns derivatives except dt.
  FullState update(const Vector3& b_acceleration, const Vector3& b_alpha, const double dt, OptionalJacobian<12, 12> F, OptionalJacobian<12, 3> G1,
                   OptionalJacobian<12, 3> G2) const;

  /// Compute tangent space contribution due to Coriolis forces
  Vector12 coriolis(double dt, const Vector3& omega, bool secondOrder = false, OptionalJacobian<12, 12> H = boost::none) const;

  /// Correct preintegrated tangent vector with our velocity and rotated gravity,
  /// taking into account Coriolis forces if omegaCoriolis is given.
  Vector12 correctPIM(const Vector12& pim, double dt, const Vector3& n_gravity, const boost::optional<Vector3>& omegaCoriolis, bool use2ndOrderCoriolis = false,
                      OptionalJacobian<12, 12> H1 = boost::none, OptionalJacobian<12, 12> H2 = boost::none) const;

  /// @}

private:
  /// @{
  /// serialization
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_NVP(R_);
    ar& BOOST_SERIALIZATION_NVP(t_);
    ar& BOOST_SERIALIZATION_NVP(v_);
    ar& BOOST_SERIALIZATION_NVP(w_);
  }
  /// @}
};

// Specialize FullState traits to use a Retract/Local that agrees with MasFactors
template <>
struct traits<FullState> : internal::Manifold<FullState>
{};

}  // namespace gtsam
