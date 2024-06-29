#include "mas_factor/mas_factor.h"

/* External or standard includes */
#include <ostream>

namespace gtsam
{

using namespace std;

//------------------------------------------------------------------------------
// Inner class PreintegratedMasMeasurements
//------------------------------------------------------------------------------

/*//{ PreintegratedMasMeasurements method */
/*//{ print() */
void PreintegratedMasMeasurements::print(const string& s) const {
  PreintegrationType::print(s);
  cout << "    preintMeasCov \n[" << preintMeasCov_ << "]" << endl;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ equals() */
bool PreintegratedMasMeasurements::equals(const PreintegratedMasMeasurements& other, double tol) const {
  return PreintegrationType::equals(other, tol) && equal_with_abs_tol(preintMeasCov_, other.preintMeasCov_, tol);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ resetIntegration() */
void PreintegratedMasMeasurements::resetIntegration() {
  PreintegrationType::resetIntegration();
  preintMeasCov_.setZero();
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ integrateMeasurement() */
void PreintegratedMasMeasurements::integrateMeasurement(const Vector& measuredMas, double dt) {

  Vector3 measuredForce(0, 0, 0);
  Vector3 measuredAcc(0, 0, 0);
  Vector3 measuredAlpha(0, 0, 0);
  Vector3 measuredTorqueDrag(0, 0, 0);
  Vector3 measuredTorqueThrust(0, 0, 0);

  MasParams p = *boost::static_pointer_cast<MasParams>(p_);

  std::vector<Vector3> rotor_p = p.getRotorPositions();

  for (int i = 0; i < measuredMas.size(); i++) {
    Vector3 motor_thrust = p.motorConstant * std::pow(measuredMas(i), 2) * Vector3(0, 0, 1);  // body frame
    measuredForce += motor_thrust;
    measuredTorqueDrag += p.rotorDirs[i] * p.momentConstant * motor_thrust;  // body frame, causes yaw
    measuredTorqueThrust += p.torqueThrustConstant * rotor_p[i].cross(motor_thrust);                    // causes roll/pitch
  }

  measuredAcc   = measuredForce / p.mass;
  measuredAlpha = p.getInertialMatrix().inverse() * (-deltaWij().cross(p.getInertialMatrix() * deltaWij()) + measuredTorqueDrag + measuredTorqueThrust);
  /* measuredAlpha = p.getInertialMatrix().inverse() * (measuredTorqueDrag + measuredTorqueThrust); */

  if (is_before_takeoff_) {
    std:: cout << "measuredMas: " << std::endl << measuredMas << std::endl;
    std:: cout << "measured_alpha: " << std::endl << measuredAlpha << std::endl;
    std:: cout << "measured_acc: " << std::endl << measuredAcc << std::endl;
    const double a_norm = measuredAcc.norm();
    const double g_norm = p.getGravity().norm();

    // TODO(petrlmat): this hack prevents the estimate to fall in Z while on the ground but should probably be handled better
    if (a_norm < g_norm) {
      measuredAcc += gtsam::Vector3(0,0,g_norm - a_norm);
    } else {
      is_before_takeoff_ = false;
    }
  }

  if (measuredAlpha.norm() > 10) {
    std:: cout << "measuredMas: " << std::endl << measuredMas << std::endl;
    std:: cout << "measured_alpha: " << std::endl << measuredAlpha << std::endl;
    std:: cout << "last_alpha: " << std::endl << last_alpha_ << std::endl;
    std:: cout << "measured_acc: " << std::endl << measuredAcc << std::endl;
    std:: cout << "last_acc_: " << std::endl << last_acc_ << std::endl;
  }

  last_acc_      = measuredAcc;
  last_acc_norm_ = measuredAcc.norm();
  last_alpha_    = measuredAlpha;
  integrateMeasurement(measuredAcc, measuredAlpha, dt);

}

void PreintegratedMasMeasurements::integrateMeasurement(const Vector3& measuredAcc, const Vector3& measuredAlpha, double dt) {
  if (dt <= 0) {
    throw std::runtime_error("PreintegratedMasMeasurements::integrateMeasurement: dt <=0");
  }

  boost::shared_ptr<MasParams> p = boost::static_pointer_cast<MasParams>(p_);

  // Update preintegrated measurements (also get Jacobian)
  Matrix_12x12 A;  // overall Jacobian wrt preintegrated measurements (df/dx)
  Matrix_12x3  B, C;
  PreintegrationType::update(measuredAcc, measuredAlpha, dt, &A, &B, &C);

  // first order covariance propagation:
  // as in [2] we consider a first order propagation that can be seen as a
  // prediction phase in EKF

  // propagate uncertainty
  const Matrix3& aCov     = p->getAccelerationCovariance();
  const Matrix3& alphaCov = p->getAlphaCovariance();
  const Matrix3& iCov     = p->getIntegrationCovariance();

  // (1/dt) allows to pass from continuous time noise to discrete time noise
  preintMeasCov_ = A * preintMeasCov_ * A.transpose();
  preintMeasCov_.noalias() += B * (aCov / dt) * B.transpose();
  preintMeasCov_.noalias() += C * (alphaCov / dt) * C.transpose();

  // NOTE(frank): (Gi*dt)*(C/dt)*(Gi'*dt), with Gi << Z_3x3, I_3x3, Z_3x3
  preintMeasCov_.block<3, 3>(0, 0).noalias() += iCov * dt;
  preintMeasCov_.block<3, 3>(3, 3).noalias() += iCov * dt;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ integrateMeasurements() */
void PreintegratedMasMeasurements::integrateMeasurements(const Matrix& measuredAccs, const Matrix& measuredAlphas, const Matrix& dts) {
  assert(measuredAccs.rows() == 3 && measuredAlphas.rows() == 3 && dts.rows() == 1);
  assert(dts.cols() >= 1);
  assert(measuredAccs.cols() == dts.cols());
  assert(measuredAlphas.cols() == dts.cols());
  size_t n = static_cast<size_t>(dts.cols());
  for (size_t j = 0; j < n; j++) {
    integrateMeasurement(measuredAccs.col(j), measuredAlphas.col(j), dts(0, j));
  }
}
/*//}*/

//------------------------------------------------------------------------------

#ifdef GTSAM_TANGENT_PREINTEGRATION
/*//{ mergeWith() */
void PreintegratedMasMeasurements::mergeWith(const PreintegratedMasMeasurements& pim12,  //
                                                    Matrix_12x12* H1, Matrix_12x12* H2) {
  PreintegrationType::mergeWith(pim12, H1, H2);
  // NOTE(gareth): Temporary P is needed as of Eigen 3.3
  const Matrix_12x12 P = *H1 * preintMeasCov_ * H1->transpose();
  preintMeasCov_       = P + *H2 * pim12.preintMeasCov_ * H2->transpose();
}
/*//}*/
#endif

//------------------------------------------------------------------------------
/*//{ deprecated */
#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
PreintegratedMasMeasurements::PreintegratedMasMeasurements(const mas_bias::ConstantBias& biasHat, const Matrix3& measuredAccCovariance,
                                                                         const Matrix3& measuredOmegaCovariance, const Matrix3& integrationErrorCovariance,
                                                                         const double mass, const double motorConstant, const double momentConstant, const double torqueThrustConstant,
                                                                         const int numRotors, const double bodyRadius, const double bodyHeight,
                                                                         const std::vector<int> rotorDirs, const bool use2ndOrderIntegration) {
  if (!use2ndOrderIntegration)
    throw("PreintegratedMasMeasurements no longer supports first-order integration: it incorrectly compensated for gravity");
  biasHat_                              = biasHat;
  boost::shared_ptr<MasParams> p = MasParams::MakeSharedU();  // TODO(petrlmat): MakeSharedD for ENU?
  p->gyroscopeCovariance                = measuredOmegaCovariance;
  p->accelerometerCovariance            = measuredAccCovariance;
  p->integrationCovariance              = integrationErrorCovariance;
  p->mass                               = mass;
  p->motorConstant                     = motorConstant;
  p->momentConstant                    = momentConstant;
  p->torqueThrustConstant              = torqueThrustConstant;
  p->numRotors                         = numRotors;
  p->bodyRadius                        = bodyRadius;
  p->bodyHeight                        = bodyHeight;
  p->rotorDirs                         = rotorDirs;
  p_                                    = p;
  resetIntegration();
}

void PreintegratedMasMeasurements::integrateMeasurement(const Vector3& measuredAcc, const Vector3& measuredAlpha, double deltaT,
                                                               boost::optional<Pose3> body_P_sensor) {
  // modify parameters to accommodate deprecated method:-(
  p_->body_P_sensor = body_P_sensor;
  integrateMeasurement(measuredAcc, measuredAlpha, deltaT);
}
#endif
/*//}*/
/*//}*/

//------------------------------------------------------------------------------
// MasFactor methods
//------------------------------------------------------------------------------

/*//{ MasFactor methods */
/*//{ MasFactor() */
MasFactor::MasFactor(Key pose_i, Key lin_vel_i, Key ang_vel_i, Key pose_j, Key lin_vel_j, Key ang_vel_j, Key bias,
                                   const PreintegratedMasMeasurements& pim)
    : Base(noiseModel::Gaussian::Covariance(pim.preintMeasCov_), pose_i, lin_vel_i, ang_vel_i, pose_j, lin_vel_j, ang_vel_j, bias), _PIM_(pim) {
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ clone() */
NonlinearFactor::shared_ptr MasFactor::clone() const {
  return boost::static_pointer_cast<NonlinearFactor>(NonlinearFactor::shared_ptr(new This(*this)));
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ operator<<() */
std::ostream& operator<<(std::ostream& os, const MasFactor& f) {
  f._PIM_.print("preintegrated measurements:\n");
  os << "  noise model sigmas: " << f.noiseModel_->sigmas().transpose();
  return os;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ print() */
void MasFactor::print(const string& s, const KeyFormatter& keyFormatter) const {
  cout << (s.empty() ? s : s + "\n") << "MasFactor(" << keyFormatter(this->key<1>())
       << "," << keyFormatter(this->key<2>()) << "," << keyFormatter(this->key<3>())
       << "," << keyFormatter(this->key<4>()) << "," << keyFormatter(this->key<5>())
       << "," << keyFormatter(this->key<6>()) << "," << keyFormatter(this->key<7>())
       << ")\n";
  cout << *this << endl;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ equals() */
bool MasFactor::equals(const NonlinearFactor& other, double tol) const {
  const This* e    = dynamic_cast<const This*>(&other);
  const bool  base = Base::equals(*e, tol);
  const bool  pim  = _PIM_.equals(e->_PIM_, tol);
  return e != nullptr && base && pim;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ evaluteError() */
Vector MasFactor::evaluateError(const Pose3& pose_i, const Vector3& lin_vel_i, const Vector3& ang_vel_i, const Pose3& pose_j, const Vector3& lin_vel_j,
                                       const Vector3& ang_vel_j, const mas_bias::ConstantBias& bias_i, boost::optional<Matrix&> H1,
                                       boost::optional<Matrix&> H2, boost::optional<Matrix&> H3, boost::optional<Matrix&> H4, boost::optional<Matrix&> H5,
                                       boost::optional<Matrix&> H6, boost::optional<Matrix&> H7) const {
  return _PIM_.computeErrorAndJacobians(pose_i, lin_vel_i, ang_vel_i, pose_j, lin_vel_j, ang_vel_j, bias_i, H1, H2, H3, H4, H5, H6, H7);
}
/*//}*/

//------------------------------------------------------------------------------

#ifdef GTSAM_TANGENT_PREINTEGRATION
/*//{ Merge() */
PreintegratedMasMeasurements MasFactor::Merge(const PreintegratedMasMeasurements& pim01,
                                                            const PreintegratedMasMeasurements& pim12) {
  if (!pim01.matchesParamsWith(pim12))
    throw std::domain_error("Cannot merge PreintegratedMasMeasurements with different params");

  if (pim01.params()->body_P_sensor)
    throw std::domain_error("Cannot merge PreintegratedMasMeasurements with sensor pose yet");

  // the bias for the merged factor will be the bias from 01
  PreintegratedMasMeasurements pim02 = pim01;

  Matrix_12x12 H1, H2;
  pim02.mergeWith(pim12, &H1, &H2);

  return pim02;
}

//------------------------------------------------------------------------------

MasFactor::shared_ptr MasFactor::Merge(const shared_ptr& f01, const shared_ptr& f12) {

  // IMU bias keys must be the same.
  if (f01->key<7>() != f12->key<7>())
    throw std::domain_error("MasFactor::Merge: IMU bias keys must be the same");

  // expect intermediate pose, velocity keys to matchup.
  if (f01->key<4>() != f12->key<1>() || f01->key<5>() != f12->key<2>() || f01->key<6>() != f12->key<3>())
    throw std::domain_error("MasFactor::Merge: intermediate pose, velocity keys need to match up");

  // return new factor
  auto pim02 = Merge(f01->preintegratedMeasurements(), f12->preintegratedMeasurements());
  return boost::make_shared<MasFactor>(f01->key<1>(),  // P0
                                              f01->key<2>(),  // V0
                                              f01->key<3>(),  // Omega0
                                              f12->key<4>(),  // P2
                                              f12->key<5>(),  // V2
                                              f12->key<6>(),  // Omega2
                                              f01->key<7>(),  // bias
                                              pim02);
}

/*//}*/
#endif

//------------------------------------------------------------------------------

/*//{ deprecated */
#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
MasFactor::MasFactor(Key pose_i, Key lin_vel_i, Key ang_vel_i, Key pose_j, Key lin_vel_j, Key ang_vel_j, Key bias,
                                   const PreintegratedMasMeasurements& pim, const Vector3& n_gravity, const Vector3& omegaCoriolis, const double mass,
                                   const double motorConstant, const double momentConstant, const double torqueThrustConstant, const int numRotors, const double bodyRadius,
                                   const double bodyHeight, const std::vector<int> rotorDirs, const boost::optional<Pose3>& body_P_sensor,
                                   const bool use2ndOrderCoriolis)
    : Base(noiseModel::Gaussian::Covariance(pim.preintMeasCov_), pose_i, lin_vel_i, ang_vel_i, pose_j, lin_vel_j, ang_vel_j, bias), _PIM_(pim) {
  boost::shared_ptr<MasParams> p = boost::make_shared<MasParams>(n_gravity);
  p->n_gravity                          = n_gravity;
  p->omegaCoriolis                      = omegaCoriolis;
  p->body_P_sensor                      = body_P_sensor;
  p->use2ndOrderCoriolis                = use2ndOrderCoriolis;
  p->mass                               = mass;
  p->motorConstant                     = motorConstant;
  p->momentConstant                    = momentConstant;
  p->torqueThrustConstant              = torqueThrustConstant;
  p->numRotors                         = numRotors;
  p->bodyRadius                        = bodyRadius;
  p->bodyHeight                        = bodyHeight;
  p->rotorDirs                         = rotorDirs;
  _PIM_.p_                              = p;
}

void MasFactor::Predict(const Pose3& pose_i, const Vector3& lin_vel_i, const Vector3& ang_vel_i, Pose3& pose_j, Vector3& lin_vel_j, Vector3& ang_vel_j,
                               const mas_bias::ConstantBias& bias_i, PreintegratedMasMeasurements& pim, const Vector3& n_gravity,
                               const Vector3& omegaCoriolis, const bool use2ndOrderCoriolis) {
  // use deprecated predict
  PoseVelocityBias pvb = pim.predict(pose_i, lin_vel_i, ang_vel_i, bias_i, n_gravity, omegaCoriolis, use2ndOrderCoriolis);
  pose_j               = pvb.pose;
  lin_vel_j            = pvb.lin_velocity;
  ang_vel_j            = pvb.ang_velocity;
}

#endif

/*//}*/ /*//}*/

//------------------------------------------------------------------------------
// MasFactor2 methods
//------------------------------------------------------------------------------

/*//{ MasFactor2 methods */
MasFactor2::MasFactor2(Key state_i, Key state_j, Key bias, const PreintegratedMasMeasurements& pim)
    : Base(noiseModel::Gaussian::Covariance(pim.preintMeasCov_), state_i, state_j, bias), _PIM_(pim) {
}

//------------------------------------------------------------------------------
NonlinearFactor::shared_ptr MasFactor2::clone() const {
  return boost::static_pointer_cast<NonlinearFactor>(NonlinearFactor::shared_ptr(new This(*this)));
}

//------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const MasFactor2& f) {
  f._PIM_.print("preintegrated measurements:\n");
  os << "  noise model sigmas: " << f.noiseModel_->sigmas().transpose();
  return os;
}

//------------------------------------------------------------------------------
void MasFactor2::print(const string& s, const KeyFormatter& keyFormatter) const {
  cout << s << "MasFactor2(" << keyFormatter(this->key1()) << "," << keyFormatter(this->key2()) << "," << keyFormatter(this->key3()) << ")\n";
  cout << *this << endl;
}

//------------------------------------------------------------------------------
bool MasFactor2::equals(const NonlinearFactor& other, double tol) const {
  const This* e    = dynamic_cast<const This*>(&other);
  const bool  base = Base::equals(*e, tol);
  const bool  pim  = _PIM_.equals(e->_PIM_, tol);
  return e != nullptr && base && pim;
}

//------------------------------------------------------------------------------
Vector MasFactor2::evaluateError(const FullState& state_i, const FullState& state_j, const mas_bias::ConstantBias& bias_i,
                                        boost::optional<Matrix&> H1, boost::optional<Matrix&> H2, boost::optional<Matrix&> H3) const {
  return _PIM_.computeError(state_i, state_j, bias_i, H1, H2, H3);
}

//------------------------------------------------------------------------------
/*//}*/

}  // namespace gtsam
