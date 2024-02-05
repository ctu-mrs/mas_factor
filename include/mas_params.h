#pragma once

#include <gtsam/navigation/PreintegrationParams.h>
#include <boost/make_shared.hpp>

namespace gtsam {

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a shared pointer to the constructor
struct GTSAM_EXPORT MasParams: PreintegrationParams {
  double mass; // total mass of UAV in kg
  double motorConstant; // the mapping of motor angular speed to thrust constant
  double momentConstant; // the mapping of motor angular speed to torque drag constant
  int numRotors; // the total number of UAV rotors
  double bodyRadius; // the radius of UAV body (without propellers) 
  double bodyHeight; // the height of UAV body (without propellers) 
  double rotorXYOffset; // the offset of rotor wrpt center of the UAV
  double rotorZOffset; // the offset of rotor wrpt center of the UAV
  std::vector<int> rotorDirs;

  bool useBiasAccCorrection = true; 
  bool useBiasAlphaCorrection = true; 

  Matrix3 accelerationCovariance;
  Matrix3 alphaCovariance;

  bool inertialMatrixComputed = false; 
  Matrix3 inertial_matrix;

  /// Default constructor for serialization only
  MasParams()
      : PreintegrationParams(),
        mass(1),
        motorConstant(10),
        momentConstant(1),
        numRotors(4),
        bodyRadius(0.25),
        bodyHeight(0.05),
        rotorXYOffset(0.1812),
        rotorZOffset(0.057),
        rotorDirs(std::vector<int>{1,1,-1,-1}) {}

  /// The Params constructor insists on getting the navigation frame gravity vector
  /// For convenience, two commonly used conventions are provided by named constructors below
  MasParams(const Vector3& n_gravity)
      : PreintegrationParams(n_gravity),
        mass(1),
        motorConstant(10),
        momentConstant(1),
        numRotors(4),
        bodyRadius(0.25),
        bodyHeight(0.05),
        rotorXYOffset(0.1812),
        rotorZOffset(0.057),
        rotorDirs(std::vector<int>{1,1,-1,-1}) {}

  // Default Params for a Z-down navigation frame, such as NED: gravity points along positive Z-axis
  static boost::shared_ptr<MasParams> MakeSharedD(double g = 9.81) {
    return boost::shared_ptr<MasParams>(new MasParams(Vector3(0, 0, g)));
  }

  // Default Params for a Z-up navigation frame, such as ENU: gravity points along negative Z-axis
  static boost::shared_ptr<MasParams> MakeSharedU(double g = 9.81) {
    return boost::shared_ptr<MasParams>(new MasParams(Vector3(0, 0, -g)));
  }

  /* void print(const std::string& s="") const; */
  /* bool equals(const PreintegratedRotationParams& other, double tol) const; */

  void setAccelerationCovariance(const Matrix3& cov) { accelerationCovariance = cov; }
  void setAlphaCovariance(const Matrix3& cov) { alphaCovariance = cov; }
  void setIntegrationCovariance(const Matrix3& cov)   { integrationCovariance = cov; }
  void setUse2ndOrderCoriolis(bool flag)              { use2ndOrderCoriolis = flag; }
  void setUseBiasAccCorrection(bool flag)              { useBiasAccCorrection = flag; }
  void setUseBiasAlphaCorrection(bool flag)              { useBiasAlphaCorrection = flag; }
  void setMass(const double m) { mass = m; inertialMatrixComputed = false;}
  void setMotorConstant(const double c_f) { motorConstant = c_f; }
  void setMomentConstant(const double c_d) { momentConstant = c_d; }
  void setNumRotors(const int n) { numRotors = n; }
  void setBodyRadius(const double r) { bodyRadius = r; inertialMatrixComputed = false;}
  void setBodyHeight(const double h) { bodyHeight = h; inertialMatrixComputed = false;}
  void setRotorDirs(const std::vector<int>  dirs) { rotorDirs = dirs; }

  const Matrix3& getAccelerationCovariance() const { return accelerationCovariance; }
  const Matrix3& getAlphaCovariance() const { return alphaCovariance; }
  const Matrix3& getIntegrationCovariance()   const { return integrationCovariance; }
  const Vector3& getGravity()   const { return n_gravity; }
  bool           getUse2ndOrderCoriolis()     const { return use2ndOrderCoriolis; }
  bool           getUseBiasAccCorrection()     const { return useBiasAccCorrection; }
  bool           getUseBiasAlphaCorrection()     const { return useBiasAlphaCorrection; }
  double getMass() const { return mass; }
  double getMotorConstant() const { return motorConstant; }
  double getMomentConstant() const { return momentConstant; }
  double getNumRotors() const { return numRotors; }
  double getFrameLength() const { return bodyRadius/sqrt(2); }
  std::vector<int> getRotorDirs() const { return rotorDirs; }
  Matrix3 getInertialMatrix() { if (!inertialMatrixComputed) { computeInertialMatrix();} return inertial_matrix;}
  std::vector<Vector3> getRotorPositions() const { return std::vector<Vector3>{Vector3(rotorXYOffset, -rotorXYOffset, rotorZOffset), Vector3(-rotorXYOffset, rotorXYOffset, rotorZOffset), Vector3(rotorXYOffset, rotorXYOffset, rotorZOffset), Vector3(-rotorXYOffset, -rotorXYOffset, rotorZOffset)}; }

  void computeInertialMatrix() {
    inertial_matrix = Z_3x3;
    inertial_matrix(0,0) = mass * (3 * pow(bodyRadius,2) + pow(bodyHeight,2)) / 12;
    inertial_matrix(1,1) = mass * (3 * pow(bodyRadius,2) + pow(bodyHeight,2)) / 12;
    inertial_matrix(2,2) = mass * pow(bodyRadius,2) / 2;
    inertialMatrixComputed = true;
  }

protected:

  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    namespace bs = ::boost::serialization;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PreintegratedRotationParams);
    ar & BOOST_SERIALIZATION_NVP(accelerationCovariance);
    ar & BOOST_SERIALIZATION_NVP(alphaCovariance);
    ar & BOOST_SERIALIZATION_NVP(integrationCovariance);
    ar & BOOST_SERIALIZATION_NVP(use2ndOrderCoriolis);
    ar & BOOST_SERIALIZATION_NVP(n_gravity);
    ar & BOOST_SERIALIZATION_NVP(mass);
    ar & BOOST_SERIALIZATION_NVP(motorConstant);
    ar & BOOST_SERIALIZATION_NVP(momentConstant);
    ar & BOOST_SERIALIZATION_NVP(numRotors);
    ar & BOOST_SERIALIZATION_NVP(bodyRadius);
    ar & BOOST_SERIALIZATION_NVP(bodyHeight);
    ar & BOOST_SERIALIZATION_NVP(inertial_matrix);
    ar & BOOST_SERIALIZATION_NVP(rotorXYOffset);
    ar & BOOST_SERIALIZATION_NVP(rotorZOffset);
    ar & BOOST_SERIALIZATION_NVP(useBiasAccCorrection);
    ar & BOOST_SERIALIZATION_NVP(useBiasAlphaCorrection);
    ar & BOOST_SERIALIZATION_NVP(rotorDirs);

  }

#ifdef GTSAM_USE_QUATERNIONS
  // Align if we are using Quaternions
public:
	GTSAM_MAKE_ALIGNED_OPERATOR_NEW
#endif
};

} // namespace gtsam
