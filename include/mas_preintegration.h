#pragma once

#include "mas_preintegration_base.h"
#include "mas_params.h"

namespace gtsam
{

/**
 * Integrate on the 12D tangent space of the FullState manifold.
 */
class GTSAM_EXPORT MasPreintegration : public MasPreintegrationBase {
protected:
  /**
   * Preintegrated navigation state, as a 12D vector on tangent space at frame i
   * Order is: theta, position, velocity
   */
  Vector12 preintegrated_;
  Vector3  last_acc_;
  double   last_acc_norm_;
  Vector3  last_alpha_;

  Eigen::Matrix<double, 12, 3> preintegrated_H_biasAcc_;
  Eigen::Matrix<double, 12, 3> preintegrated_H_biasAlpha_;

  /// Default constructor for serialization
  MasPreintegration() {
    resetIntegration();
  }

public:
  /// @name Constructors/destructors
  /// @{

  /**
   *  Constructor, initializes the variables in the base class
   *  @param p    Parameters, typically fixed in a single application
   *  @param bias Current estimate of acceleration and rotation rate biases
   */
  MasPreintegration(const boost::shared_ptr<Params>& p, const mas_bias::ConstantBias& biasHat = mas_bias::ConstantBias());

  /// Virtual destructor
  virtual ~MasPreintegration() {
  }

  /// @}

  /// @name Basic utilities
  /// @{
  /// Re-initialize PreintegratedMeasurements
  void resetIntegration() override;

  /// @}

  /// @name Instance variables access
  /// @{
  Vector3 deltaPij() const override {
    return preintegrated_.segment<3>(3);
  }
  Vector3 deltaVij() const override {
    return preintegrated_.segment<3>(6);
  }
  Vector3 deltaWij() const override {
    return preintegrated_.segment<3>(9);
  }
  Rot3 deltaRij() const override {
    return Rot3::Expmap(theta());
  }
  FullState deltaXij() const override {
    return FullState().retract(preintegrated_);
  }
  Vector3 lastAcc() const {
    return last_acc_;
  }
  double lastAccNorm() const {
    return last_acc_norm_;
  }
  Vector3 lastAlpha() const {
    return last_alpha_;
  }

  const Vector12& preintegrated() const {
    return preintegrated_;
  }
  Vector3 theta() const {
    return preintegrated_.head<3>();
  }
  const Eigen::Matrix<double, 12, 3>& preintegrated_H_biasAcc() const {
    return preintegrated_H_biasAcc_;
  }
  const Eigen::Matrix<double, 12, 3>& preintegrated_H_biasAlpha() const {
    return preintegrated_H_biasAlpha_;
  }


  const Vector12& getPreintegrated() const override {
    return preintegrated();
  }

  void setAccelerationCovariance(const gtsam::Matrix33& cov) {
    boost::static_pointer_cast<MasParams>(p_)->setAccelerationCovariance(cov) ;
  }

  /// @name Testable
  /// @{
  bool equals(const MasPreintegration& other, double tol) const;
  /// @}

  /// @name Main functionality
  /// @{

  // Update integrated vector on tangent manifold preintegrated with acceleration
  // Static, functional version.
  static Vector12 UpdatePreintegrated(const Vector3& a_body, const Vector3& alpha_body, const double dt, const Vector12& preintegrated,
                                      OptionalJacobian<12, 12> A = boost::none, OptionalJacobian<12, 3> B = boost::none,
                                      OptionalJacobian<12, 3> C = boost::none);

  /// Update preintegrated measurements and get derivatives
  /// It takes measured quantities in the j frame
  /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
  void update(const Vector3& measuredAcc, const Vector3& measuredAlpha, const double dt, Eigen::Matrix<double, 12, 12>* A, Eigen::Matrix<double, 12, 3>* B,
              Eigen::Matrix<double, 12, 3>* C);

  /// Given the estimate of the bias, return a FullState tangent vector
  /// summarizing the preintegrated MAS measurements so far
  Vector12 biasCorrectedDelta(const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 6> H = boost::none) const;


  // Compose the two pre-integrated 9D-vectors zeta01 and zeta02, with derivatives
  static Vector12 Compose(const Vector12& zeta01, const Vector12& zeta12, double deltaT12, OptionalJacobian<12, 12> H1 = boost::none,
                          OptionalJacobian<12, 12> H2 = boost::none);

  /// Merge in a different set of measurements and update bias derivatives accordingly
  /// The derivatives apply to the preintegrated Vector12
  void mergeWith(const MasPreintegration& pim, Matrix_12x12* H1, Matrix_12x12* H2);
  /// @}

  /** Dummy clone for MATLAB */
  virtual boost::shared_ptr<MasPreintegration> clone() const {
    return boost::shared_ptr<MasPreintegration>();
  }

  /// @}

private:
  /** Serialization function */
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    namespace bs = ::boost::serialization;
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(MasPreintegrationBase);
    ar& BOOST_SERIALIZATION_NVP(preintegrated_);
    ar & BOOST_SERIALIZATION_NVP(preintegrated_H_biasAcc_);
    ar & BOOST_SERIALIZATION_NVP(preintegrated_H_biasAlpha_);
  }

public:
  GTSAM_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace gtsam
