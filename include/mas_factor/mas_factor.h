#pragma once

/* GTSAM includes */
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/base/debug.h>

#include "mas_params.h"
#include "mas_preintegration.h"

namespace gtsam {

/*//{ class NoiseModelFactor7*/
template<class VALUE1, class VALUE2, class VALUE3, class VALUE4, class VALUE5, class VALUE6, class VALUE7>
class NoiseModelFactor7: public NoiseModelFactor {

public:

  // typedefs for value types pulled from keys
  typedef VALUE1 X1;
  typedef VALUE2 X2;
  typedef VALUE3 X3;
  typedef VALUE4 X4;
  typedef VALUE5 X5;
  typedef VALUE6 X6;
  typedef VALUE7 X7;

protected:

  typedef NoiseModelFactor Base;
  typedef NoiseModelFactor7<VALUE1, VALUE2, VALUE3, VALUE4, VALUE5, VALUE6, VALUE7> This;

public:

  /**
   * Default Constructor for I/O
   */
  NoiseModelFactor7() {}

  /**
   * Constructor
   * @param noiseModel shared pointer to noise model
   * @param j1 key of the first variable
   * @param j2 key of the second variable
   * @param j3 key of the third variable
   * @param j4 key of the fourth variable
   * @param j5 key of the fifth variable
   * @param j6 key of the fifth variable
   */
  NoiseModelFactor7(const SharedNoiseModel& noiseModel, Key j1, Key j2, Key j3, Key j4, Key j5, Key j6, Key j7) :
    Base(noiseModel, cref_list_of<7>(j1)(j2)(j3)(j4)(j5)(j6)(j7)) {}

  virtual ~NoiseModelFactor7() {}

  /** methods to retrieve keys */
  inline Key key1() const { return keys_[0]; }
  inline Key key2() const { return keys_[1]; }
  inline Key key3() const { return keys_[2]; }
  inline Key key4() const { return keys_[3]; }
  inline Key key5() const { return keys_[4]; }
  inline Key key6() const { return keys_[5]; }
  inline Key key7() const { return keys_[6]; }

  /** Calls the 6-key specific version of evaluateError, which is pure virtual
   * so must be implemented in the derived class. */
  virtual Vector unwhitenedError(const Values& x, boost::optional<std::vector<Matrix>&> H = boost::none) const {
    if(this->active(x)) {
      if(H)
        return evaluateError(x.at<X1>(keys_[0]), x.at<X2>(keys_[1]), x.at<X3>(keys_[2]), x.at<X4>(keys_[3]), x.at<X5>(keys_[4]), x.at<X6>(keys_[5]), x.at<X7>(keys_[6]), (*H)[0], (*H)[1], (*H)[2], (*H)[3], (*H)[4], (*H)[5], (*H)[6]);
      else
        return evaluateError(x.at<X1>(keys_[0]), x.at<X2>(keys_[1]), x.at<X3>(keys_[2]), x.at<X4>(keys_[3]), x.at<X5>(keys_[4]), x.at<X6>(keys_[5]), x.at<X7>(keys_[6]));
    } else {
      return Vector::Zero(this->dim());
    }
  }

  /**
   *  Override this method to finish implementing a 6-way factor.
   *  If any of the optional Matrix reference arguments are specified, it should compute
   *  both the function evaluation and its derivative(s) in X1 (and/or X2, X3).
   */
  virtual Vector
  evaluateError(const X1&, const X2&, const X3&, const X4&, const X5&, const X6&, const X7&,
      boost::optional<Matrix&> H1 = boost::none,
      boost::optional<Matrix&> H2 = boost::none,
      boost::optional<Matrix&> H3 = boost::none,
      boost::optional<Matrix&> H4 = boost::none,
      boost::optional<Matrix&> H5 = boost::none,
      boost::optional<Matrix&> H6 = boost::none,
      boost::optional<Matrix&> H7 = boost::none) const = 0;

private:

  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar & boost::serialization::make_nvp("NoiseModelFactor",
        boost::serialization::base_object<Base>(*this));
  }
};
/*//}*/

/*//{ class PreintegratedMasMeasurements */

typedef MasPreintegration PreintegrationType;

class GTSAM_EXPORT PreintegratedMasMeasurements: public PreintegrationType {

  friend class MasFactor;
  friend class MasFactor2;

protected:

  Matrix_12x12 preintMeasCov_; ///< COVARIANCE OF: [PreintROTATION PreintPOSITION PreintVELOCITY PreintOMEGA]
  ///< (first-order propagation from *measurementCovariance*).

  bool is_before_takeoff_ = true;

public:

  /// Default constructor for serialization and Cython wrapper
  PreintegratedMasMeasurements() {
    preintMeasCov_.setZero();
  }

 /**
   *  Constructor, initializes the class with no measurements
   *  @param p       Parameters, typically fixed in a single application
   *  @param biasHat Current estimate of acceleration and rotation rate biases
   */
  PreintegratedMasMeasurements(const boost::shared_ptr<MasParams>& p, const mas_bias::ConstantBias& biasHat = mas_bias::ConstantBias()) :
      PreintegrationType(p) {
    preintMeasCov_.setZero();
  }

/**
  *  Construct preintegrated directly from members: base class and preintMeasCov
  *  @param base               PreintegrationType instance
  *  @param preintMeasCov      Covariance matrix used in noise model.
  */
  PreintegratedMasMeasurements(const PreintegrationType& base, const Matrix_12x12& preintMeasCov)
     : PreintegrationType(base),
       preintMeasCov_(preintMeasCov) {
  }

  /// Virtual destructor
  virtual ~PreintegratedMasMeasurements() {
  }

  /// print
  void print(const std::string& s = "Preintegrated Measurements:") const override;

  /// equals
  bool equals(const PreintegratedMasMeasurements& expected, double tol = 1e-9) const;

  /// Re-initialize PreintegratedMasMeasurements
  void resetIntegration() override;

  /**
   * Add a single motor angular speed measurement to the preintegration.
   * @param measuredMass Measured motor angular speeds
   * @param dt Time interval between this and the last motor angular speed measurement
   */
  void integrateMeasurement(
    const Vector3& measuredAcc, const Vector3& measuredAlpha, double dt) override;

  void integrateMeasurement(const Vector& measuredMass,
       const double dt);

  /// Add multiple measurements, in matrix columns
  void integrateMeasurements(const Matrix& measuredAccs, const Matrix& measuredAlphas,
                             const Matrix& dts);

  /// Return pre-integrated measurement covariance
  Matrix preintMeasCov() const { return preintMeasCov_; }

#ifdef GTSAM_TANGENT_PREINTEGRATION
  /// Merge in a different set of measurements and update bias derivatives accordingly
  void mergeWith(const PreintegratedMasMeasurements& pim, Matrix_12x12* H1, Matrix_12x12* H2);
#endif

#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
  /// @deprecated constructor
  /// NOTE(frank): assumes Z-Down convention, only second order integration supported
  PreintegratedMasMeasurements(
      const mas_bias::ConstantBias& biasHat,
      const Matrix3& measuredAccCovariance,
      const Matrix3& measuredAlphaCovariance,
      const Matrix3& integrationErrorCovariance,
      const double mass,
      const double motor_constant,
      const double moment_constant,
      const int num_rotors,
      const double body_radius,
      const double body_height,
      const std::vector<int> rotor_dirs,
      bool use2ndOrderIntegration = true);

  /// @deprecated version of integrateMeasurement with body_P_sensor
  /// Use parameters instead
  void integrateMeasurement(const Vector3& measuredAcc,
      const Vector3& measuredAlpha, double dt,
      boost::optional<Pose3> body_P_sensor);
#endif

private:

  /// Serialization function
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    namespace bs = ::boost::serialization;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PreintegrationType);
    ar & BOOST_SERIALIZATION_NVP(preintMeasCov_);
  }
};
/*//}*/

/*//{ class MasFactor */
class GTSAM_EXPORT MasFactor: public NoiseModelFactor7<Pose3, Vector3, Vector3, Pose3, Vector3, Vector3, mas_bias::ConstantBias> {
private:

  typedef MasFactor This;
  typedef NoiseModelFactor7<Pose3, Vector3, Vector3, Pose3, Vector3, Vector3, mas_bias::ConstantBias> Base;

  PreintegratedMasMeasurements _PIM_;

public:

  /** Shorthand for a smart pointer to a factor */
#if !defined(_MSC_VER) && __GNUC__ == 4 && __GNUC_MINOR__ > 5
  typedef typename boost::shared_ptr<MasFactor> shared_ptr;
#else
  typedef boost::shared_ptr<MasFactor> shared_ptr;
#endif

  /** Default constructor - only use for serialization */
  MasFactor() {}

  /**
   * Constructor
   * @param pose_i Previous pose key
   * @param lin_vel_i  Previous linear velocity key
   * @param ang_vel_i  Previous angular velocity key
   * @param pose_j Current pose key
   * @param lin_vel_j  Current linear velocity key
   * @param ang_vel_j  Current angular velocity key
   * @param bias   Previous bias key
   */
  MasFactor(Key pose_i, Key lin_vel_i, Key ang_vel_i, Key pose_j, Key lin_vel_j, Key ang_vel_j, Key bias,
      const PreintegratedMasMeasurements& preintegratedMeasurements);

  virtual ~MasFactor() {
  }

  /// @return a deep copy of this factor
  virtual gtsam::NonlinearFactor::shared_ptr clone() const;

  /// @name Testable
  /// @{
  GTSAM_EXPORT friend std::ostream& operator<<(std::ostream& os, const MasFactor&);
  virtual void print(const std::string& s, const KeyFormatter& keyFormatter =
      DefaultKeyFormatter) const;
  virtual bool equals(const NonlinearFactor& expected, double tol = 1e-9) const;
  /// @}

  /** Access the preintegrated measurements. */

  const PreintegratedMasMeasurements& preintegratedMeasurements() const {
    return _PIM_;
  }

  /** implement functions needed to derive from Factor */

  /// vector of errors
  Vector evaluateError(const Pose3& pose_i, const Vector3& lin_vel_i, const Vector3& ang_vel_i,
      const Pose3& pose_j, const Vector3& lin_vel_j, const Vector3& ang_vel_j, const mas_bias::ConstantBias& bias_i,
      boost::optional<Matrix&> H1 =
          boost::none, boost::optional<Matrix&> H2 = boost::none,
      boost::optional<Matrix&> H3 = boost::none, boost::optional<Matrix&> H4 =
          boost::none, boost::optional<Matrix&> H5 = boost::none, boost::optional<Matrix&> H6 = boost::none, boost::optional<Matrix&> H7 = boost::none) const;

#ifdef GTSAM_TANGENT_PREINTEGRATION
  /// Merge two pre-integrated measurement classes
  static PreintegratedMasMeasurements Merge(
      const PreintegratedMasMeasurements& pim01,
      const PreintegratedMasMeasurements& pim12);

  /// Merge two factors
  static shared_ptr Merge(const shared_ptr& f01, const shared_ptr& f12);
#endif

#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
  /// @deprecated typename
  typedef PreintegratedMasMeasurements PreintegratedMeasurements;

  /// @deprecated constructor, in the new one gravity, coriolis settings are in PreintegrationParams
  MasFactor(Key pose_i, Key lin_vel_i, Key ang_vel_i, Key pose_j, Key lin_vel_j, Key ang_vel_j, Key bias,
      const PreintegratedMeasurements& preintegratedMeasurements,
      const Vector3& n_gravity, const Vector3& omegaCoriolis, const double mass, const double motor_constant, const double moment_constant, const int num_rotors, const double body_radius, const double body_height, const std::vector<int> rotor_dirs,
      const boost::optional<Pose3>& body_P_sensor = boost::none,
      const bool use2ndOrderCoriolis = false);

  /// @deprecated use PreintegrationBase::predict,
  /// in the new one gravity, coriolis settings are in PreintegrationParams
  static void Predict(const Pose3& pose_i, const Vector3& lin_vel_i, const Vector3& ang_vel_i, Pose3& pose_j,
      Vector3& lin_vel_j, Vector3& ang_vel_j, const mas_bias::ConstantBias& bias_i,
      PreintegratedMeasurements& pim, const Vector3& n_gravity,
      const Vector3& omegaCoriolis, const bool use2ndOrderCoriolis = false);
#endif

private:

  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar & boost::serialization::make_nvp("NoiseModelFactor7",
         boost::serialization::base_object<Base>(*this));
    ar & BOOST_SERIALIZATION_NVP(_PIM_);
  }
};
// class MasFactor
/*//}*/

class GTSAM_EXPORT MasFactor2 : public NoiseModelFactor3<FullState, FullState, mas_bias::ConstantBias> {
private:

  typedef MasFactor2 This;
  typedef NoiseModelFactor3<FullState, FullState, mas_bias::ConstantBias> Base;

  PreintegratedMasMeasurements _PIM_;

public:

  /** Default constructor - only use for serialization */
  MasFactor2() {}

  /**
   * Constructor
   * @param state_i Previous state key
   * @param state_j Current state key
   * @param bias    Previous bias key
   */
  MasFactor2(Key state_i, Key state_j, Key bias,
             const PreintegratedMasMeasurements& preintegratedMeasurements);

  virtual ~MasFactor2() {
  }

  /// @return a deep copy of this factor
  virtual gtsam::NonlinearFactor::shared_ptr clone() const;

  /// @name Testable
  /// @{
  GTSAM_EXPORT friend std::ostream& operator<<(std::ostream& os, const MasFactor2&);
  virtual void print(const std::string& s, const KeyFormatter& keyFormatter =
      DefaultKeyFormatter) const;
  virtual bool equals(const NonlinearFactor& expected, double tol = 1e-9) const;
  /// @}

  /** Access the preintegrated measurements. */

  const PreintegratedMasMeasurements& preintegratedMeasurements() const {
    return _PIM_;
  }

  /** implement functions needed to derive from Factor */

  /// vector of errors
  Vector evaluateError(const FullState& state_i, const FullState& state_j, const mas_bias::ConstantBias& bias_i,
                       boost::optional<Matrix&> H1 = boost::none,
                       boost::optional<Matrix&> H2 = boost::none,
                       boost::optional<Matrix&> H3 = boost::none
                       ) const;

private:

  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar & boost::serialization::make_nvp("NoiseModelFactor3",
         boost::serialization::base_object<Base>(*this));
    ar & BOOST_SERIALIZATION_NVP(_PIM_);
  }
};
// class MasFactor2

template <>
struct traits<PreintegratedMasMeasurements> : public Testable<PreintegratedMasMeasurements> {};

template <>
struct traits<MasFactor> : public Testable<MasFactor> {};

template <>
struct traits<MasFactor2> : public Testable<MasFactor2> {};

} /// namespace gtsam