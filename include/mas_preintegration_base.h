#pragma once

#include <gtsam/navigation/PreintegrationParams.h>
#include <gtsam/linear/NoiseModel.h>

#include <iosfwd>
#include <string>
#include <utility>

#include "full_state.h"
#include "mas_factor/mas_bias.h"

namespace gtsam
{

typedef Eigen::Matrix<double, 12, 12> Matrix_12x12;
typedef Eigen::Matrix<double, 12, 6>  Matrix_12x6;
typedef Eigen::Matrix<double, 12, 3>  Matrix_12x3;

typedef Eigen::Matrix<double, 12, 12> Matrix_12x12;

 static const Eigen::MatrixBase<Matrix_12x3>::ConstantReturnType Z_12x3 = Matrix_12x3::Zero();

#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
/// @deprecated
struct PoseVelocityBias
{
  Pose3                 pose;
  Vector3               lin_velocity;
  Vector3               ang_velocity;
  mas_bias::ConstantBias bias;
  PoseVelocityBias(const Pose3& _pose, const Vector3& _lin_velocity, const Vector3& _ang_velocity, const mas_bias::ConstantBias _bias)
      : pose(_pose), lin_velocity(_lin_velocity), ang_velocity(_ang_velocity), bias(_bias) {
  }
  PoseVelocityBias(const FullState& fullState, const mas_bias::ConstantBias _bias)
      : pose(fullState.pose()), lin_velocity(fullState.linVelocity()), ang_velocity(fullState.angVelocity()), bias(_bias) {
  }
  FullState navState() const {
    return FullState(pose, lin_velocity, ang_velocity);
  }
};
#endif

/**
 * MasPreintegrationBase is the base class for PreintegratedMeasurements
 * (in MasFactor). It includes the definitions of the preintegrated variables and the methods to access, print, and compare them.
 */
class GTSAM_EXPORT MasPreintegrationBase {
public:
  typedef PreintegrationParams Params;
  typedef mas_bias::ConstantBias Bias;

protected:
  /// Parameters. Declared mutable only for deprecated predict method.
#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
  mutable
#endif
      boost::shared_ptr<Params>
          p_;

  /// Acceleration and gyro bias used for preintegration
  Bias biasHat_;

  /// Time interval from i to j
  double deltaTij_;

  /// Default constructor for serialization
  MasPreintegrationBase() {
  }

  /// Virtual destructor for serialization
  virtual ~MasPreintegrationBase() {
  }

public:
  /// @name Constructors
  /// @{

  /**
   *  Constructor, initializes the variables in the base class
   *  @param p    Parameters, typically fixed in a single application
   *  @param bias Current estimate of acceleration and rotation rate biases
   */
  MasPreintegrationBase(const boost::shared_ptr<Params>& p, const mas_bias::ConstantBias& biasHat = mas_bias::ConstantBias());

  /// @}

  /// @name Basic utilities
  /// @{
  /// Re-initialize PreintegratedMeasurements
  virtual void resetIntegration() = 0;

  /// @name Basic utilities
  /// @{
  /// Re-initialize PreintegratedMeasurements and set new bias
  void resetIntegrationAndSetBias(const Bias& biasHat);

  /// check parameters equality: checks whether shared pointer points to same Params object.
  bool matchesParamsWith(const MasPreintegrationBase& other) const {
    return p_.get() == other.p_.get();
  }

  /// shared pointer to params
  const boost::shared_ptr<Params>& params() const {
    return p_;
  }

  /// const reference to params
  const Params& p() const {
    return *boost::static_pointer_cast<Params>(p_);
  }

#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
  void set_body_P_sensor(const Pose3& body_P_sensor) {
    p_->body_P_sensor = body_P_sensor;
  }
#endif
  /// @}

  /// @name Instance variables access
  /// @{
  const mas_bias::ConstantBias& biasHat() const {
    return biasHat_;
  }
  double deltaTij() const {
    return deltaTij_;
  }

  virtual Vector3   deltaPij() const = 0;
  virtual Vector3   deltaVij() const = 0;
  virtual Vector3   deltaWij() const = 0;
  virtual Rot3      deltaRij() const = 0;
  virtual FullState deltaXij() const = 0;

  virtual const Vector12& getPreintegrated() const = 0;

  /// @}

  /// @name Testable
  /// @{
  GTSAM_EXPORT friend std::ostream& operator<<(std::ostream& os, const MasPreintegrationBase& pim);
  virtual void                      print(const std::string& s) const;
  /// @}

  /// @name Main functionality
  /// @{

  /**
   * Subtract estimate and correct for sensor pose
   * Compute the derivatives due to non-identity body_P_sensor (rotation and centrifugal acc)
   * Ignore D_correctedAlpha_measuredAcc as it is trivially zero
   */
  /* std::pair<Vector3, Vector3> correctMeasurementsBySensorPose(const Vector3& unbiasedAcc, const Vector3& unbiasedAlpha, */
  /*                                                             OptionalJacobian<3, 3> correctedAcc_H_unbiasedAcc     = boost::none, */
  /*                                                             OptionalJacobian<3, 3> correctedAcc_H_unbiasedAlpha   = boost::none, */
  /*                                                             OptionalJacobian<3, 3> correctedAlpha_H_unbiasedAlpha = boost::none) const; */

  /**
   *  Update preintegrated measurements and get derivatives
   * It takes measured quantities in the j frame
   * Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
   */
  virtual void update(const Vector3& measuredAcc, const Vector3& measuredAlpha, const double dt, Matrix_12x12* A, Matrix_12x3* B, Matrix_12x3* C) = 0;

  /// Version without derivatives
  virtual void integrateMeasurement(const Vector3& measuredAcc, const Vector3& measuredAlpha, const double dt);

  /// Given the estimate of the bias, return a FullState tangent vector
  /// summarizing the preintegrated MAS measurements so far
  virtual Vector12 biasCorrectedDelta(const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 6> H = boost::none) const = 0;

  /// Predict state at time j
  FullState predict(const FullState& state_i, const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 12> H1 = boost::none,
                    OptionalJacobian<12, 6> H2 = boost::none) const;

  /// Calculate error given fullStates
  Vector12 computeError(const FullState& state_i, const FullState& state_j, const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 12> H1,
                        OptionalJacobian<12, 12> H2, OptionalJacobian<12, 6> H3) const;


  /**
   * Compute errors w.r.t. preintegrated measurements and jacobians
-   * wrt pose_i, vel_i, bias_i, pose_j, bias_j
   */
  Vector12 computeErrorAndJacobians(const Pose3& pose_i, const Vector3& vel_i, const Vector3& omega_i, const Pose3& pose_j, const Vector3& vel_j,
                                    const Vector3& omega_j, const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 6> H1 = boost::none,
                                    OptionalJacobian<12, 3> H2 = boost::none, OptionalJacobian<12, 3> H3 = boost::none,
                                    OptionalJacobian<12, 6> H4 = boost::none, OptionalJacobian<12, 3> H5 = boost::none,
                                    OptionalJacobian<12, 3> H6 = boost::none, OptionalJacobian<12, 6> H7 = boost::none) const;

#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
  /// @name Deprecated
  /// @{

  /// @deprecated predict
  PoseVelocityBias predict(const Pose3& pose_i, const Vector3& lin_vel_i, const Vector3& ang_vel_i, const mas_bias::ConstantBias& bias_i,
                           const Vector3& n_gravity, const Vector3& omegaCoriolis, const bool use2ndOrderCoriolis = false) const;

  /// @}
#endif

private:
  /** Serialization function */
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_NVP(p_);
    ar& BOOST_SERIALIZATION_NVP(deltaTij_);
    ar& BOOST_SERIALIZATION_NVP(biasHat_);
  }

public:
  GTSAM_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace gtsam
