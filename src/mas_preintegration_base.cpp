#include "mas_preintegration_base.h"
#include <gtsam/base/numericalDerivative.h>
#include <boost/make_shared.hpp>

using namespace std;

namespace gtsam
{

/*//{ constructor */
MasPreintegrationBase::MasPreintegrationBase(const boost::shared_ptr<Params>& p, const Bias& biasHat) : p_(p), biasHat_(biasHat), deltaTij_(0.0) {
}
/*//}*/

/*//{ operator<<() */
ostream& operator<<(ostream& os, const MasPreintegrationBase& pim) {
  os << "    deltaTij " << pim.deltaTij_ << endl;
  os << "    deltaRij.ypr = (" << pim.deltaRij().ypr().transpose() << ")" << endl;
  os << "    deltaPij " << Point3(pim.deltaPij()) << endl;
  os << "    deltaVij " << Point3(pim.deltaVij()) << endl;
  os << "    deltaWij " << Point3(pim.deltaWij()) << endl;
  os << "    biasAcc " << Point3(pim.biasHat_.linAcc()) << endl;
  os << "    biasAlpha " << Point3(pim.biasHat_.angAcc()) << endl;
  return os;
}
/*//}*/

/*//{ print() */
void MasPreintegrationBase::print(const string& s) const {
  cout << s << *this << endl;
}
/*//}*/

/*//{ resetIntegrationAndSetBias() */
void MasPreintegrationBase::resetIntegrationAndSetBias(const Bias& biasHat) {
  biasHat_ = biasHat;
  resetIntegration();
}
/*//}*/

/*//{ integrateMeasurement() */
void MasPreintegrationBase::integrateMeasurement(const Vector3& measuredAcc, const Vector3& measuredAlpha, double dt) {
  Matrix_12x12 A;
  Matrix_12x3  B, C;
  update(measuredAcc, measuredAlpha, dt, &A, &B, &C);
}
/*//}*/

/*//{ predict() */
FullState MasPreintegrationBase::predict(const FullState& state_i, const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 12> H1,
                                                OptionalJacobian<12, 6> H2) const {
  Matrix_12x6 D_biasCorrected_bias;

  Vector12 biasCorrected = biasCorrectedDelta(bias_i, H2 ? &D_biasCorrected_bias : nullptr);

  // Correct for initial velocity and gravity
  Matrix_12x12 D_delta_state, D_delta_biasCorrected;
  Vector12     xi = state_i.correctPIM(biasCorrected, deltaTij_, p().n_gravity, p().omegaCoriolis, p().use2ndOrderCoriolis, H1 ? &D_delta_state : nullptr,
                                   H2 ? &D_delta_biasCorrected : nullptr);


  // Use retract to get back to FullState manifold
  Matrix_12x12 D_predict_state, D_predict_delta;

  FullState state_j = state_i.retract(xi, H1 ? &D_predict_state : nullptr, H2 || H2 ? &D_predict_delta : nullptr);
  /* cout << "state_i:" << endl; */
  /* state_i.print(); */
  /* cout << "xi:" << endl; */
  /* cout << xi.transpose() << endl; */
  /* cout << "state_j:" << endl; */
  /* state_j.print(); */

  if (H1) {
    *H1 = D_predict_state + D_predict_delta * D_delta_state;
    /* cout << "D_predict_state:" << endl << D_predict_state << endl; */
    /* cout << "+ D_predict_delta:" << endl << D_predict_delta << endl; */
    /* cout << "* D_delta_state:" << endl << D_delta_state << endl; */
    /* cout << "H:" << endl << *H1 << endl; */
  }

  if (H2) {
    *H2 = D_predict_delta * D_delta_biasCorrected * D_biasCorrected_bias;
  }

  return state_j;
}
/*//}*/

/*//{ computeError() */
Vector12 MasPreintegrationBase::computeError(const FullState& state_i, const FullState& state_j, const mas_bias::ConstantBias& bias_i,
                                                    OptionalJacobian<12, 12> H1, OptionalJacobian<12, 12> H2, OptionalJacobian<12, 6> H3) const {
  // Predict state at time j
  Matrix_12x12 D_predict_state_i;
  Matrix_12x6  D_predict_bias_i;
  FullState    predictedState_j = predict(state_i, bias_i, H1 ? &D_predict_state_i : 0, H3 ? &D_predict_bias_i : 0);

  // Calculate error
  Matrix_12x12 D_error_state_j, D_error_predict;
  Vector12     error = state_j.localCoordinates(predictedState_j, H2 ? &D_error_state_j : 0, H1 || H3 ? &D_error_predict : 0);

  if (H1)
    *H1 << D_error_predict * D_predict_state_i;
  if (H2)
    *H2 << D_error_state_j;
  if (H3)
    *H3 << D_error_predict * D_predict_bias_i;

  return error;
}
/*//}*/

/*//{ computeErrorAndJacobians() */
Vector12 MasPreintegrationBase::computeErrorAndJacobians(const Pose3& pose_i, const Vector3& vel_i, const Vector3& omega_i, const Pose3& pose_j,
                                                                const Vector3& vel_j, const Vector3& omega_j, const mas_bias::ConstantBias& bias_i,
                                                                OptionalJacobian<12, 6> H1, OptionalJacobian<12, 3> H2, OptionalJacobian<12, 3> H3,
                                                                OptionalJacobian<12, 6> H4, OptionalJacobian<12, 3> H5, OptionalJacobian<12, 3> H6,
                                                                OptionalJacobian<12, 6> H7) const {

  // Note that derivative of constructors below is not identity for velocity, but
  // a 9*3 matrix == Z_3x3, Z_3x3, state.R().transpose()
  FullState state_i(pose_i, vel_i, omega_i);
  FullState state_j(pose_j, vel_j, omega_j);

  // Predict state at time j
  Matrix_12x12 D_error_state_i, D_error_state_j;
  Vector12     error = computeError(state_i, state_j, bias_i, H1 || H2 || H3 ? &D_error_state_i : 0, H4 || H5 || H6 ? &D_error_state_j : 0, H7);

  // Separate out derivatives in terms of 5 arguments
  // Note that doing so requires special treatment of velocities, as when treated as
  // separate variables the retract applied will not be the semi-direct product in FullState
  // Instead, the velocities in nav are updated using a straight addition
  // This is difference is accounted for by the R().transpose calls below
  if (H1)
    *H1 << D_error_state_i.leftCols<6>();
  if (H2)
    *H2 << D_error_state_i.block<12, 3>(0, 6) * state_i.R().transpose();
  if (H3)
    *H3 << D_error_state_i.rightCols<3>() * state_i.R().transpose();
  if (H4)
    *H4 << D_error_state_j.leftCols<6>();
  if (H5)
    *H5 << D_error_state_j.block<12, 3>(0, 6) * state_j.R().transpose();
  if (H6)
    *H6 << D_error_state_j.rightCols<3>() * state_j.R().transpose();
  // NOTE(petrlmat): H7 was missing in bias version

  return error;
}
/*//}*/

#ifdef GTSAM_ALLOW_DEPRECATED_SINCE_V4
/*//{ predict() */
PoseVelocityBias MasPreintegrationBase::predict(const Pose3& pose_i, const Vector3& lin_vel_i, const Vector3& ang_vel_i,
                                                       const mas_bias::ConstantBias& bias_i, const Vector3& n_gravity, const Vector3& omegaCoriolis,
                                                       const bool use2ndOrderCoriolis) const {
  cout << "USING DEPRECATED FUNCTION!" << endl;
  // NOTE(frank): parameters are supposed to be constant, below is only provided for compatibility
  boost::shared_ptr<Params> q = boost::make_shared<Params>(p());
  q->n_gravity                = n_gravity;
  q->omegaCoriolis            = omegaCoriolis;
  q->use2ndOrderCoriolis      = use2ndOrderCoriolis;
  p_                          = q;
  return PoseVelocityBias(predict(FullState(pose_i, lin_vel_i, ang_vel_i), bias_i), bias_i);
}
/*//}*/
#endif

}  // namespace gtsam
