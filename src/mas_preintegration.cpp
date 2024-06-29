#include "mas_preintegration.h"
#include <gtsam/base/numericalDerivative.h>
#include <boost/make_shared.hpp>

using namespace std;

namespace gtsam
{

//------------------------------------------------------------------------------

/*//{ MasPreintegration() */
MasPreintegration::MasPreintegration(const boost::shared_ptr<Params>& p, const Bias& biasHat) : MasPreintegrationBase(p, biasHat) {
  resetIntegration();
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ resetIntegration() */
void MasPreintegration::resetIntegration() {
  deltaTij_ = 0.0;
  preintegrated_.setZero();
  preintegrated_H_biasAcc_.setZero();
  preintegrated_H_biasAlpha_.setZero();
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ equals() */
bool MasPreintegration::equals(const MasPreintegration& other, double tol) const {
  return p_->equals(*other.p_, tol) && std::abs(deltaTij_ - other.deltaTij_) < tol && biasHat_.equals(other.biasHat_, tol) &&
         equal_with_abs_tol(preintegrated_, other.preintegrated_, tol) && equal_with_abs_tol(preintegrated_H_biasAcc_, other.preintegrated_H_biasAcc_, tol) &&
         equal_with_abs_tol(preintegrated_H_biasAlpha_, other.preintegrated_H_biasAlpha_, tol);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ UpdatePreintegrated() */
Vector12 MasPreintegration::UpdatePreintegrated(const Vector3& a_body, const Vector3& alpha_body, double dt, const Vector12& preintegrated,
                                                       OptionalJacobian<12, 12> A, OptionalJacobian<12, 3> B, OptionalJacobian<12, 3> C) {

  const auto theta    = preintegrated.segment<3>(0);
  const auto position = preintegrated.segment<3>(3);
  const auto velocity = preintegrated.segment<3>(6);
  const auto omega    = preintegrated.segment<3>(9);

  // This functor allows for saving computation when exponential map and its
  // derivatives are needed at the same location in so<3>
  so3::DexpFunctor local(theta);

  // Calculate exact mean propagation
  Matrix3       alpha_tangent_H_theta, invH;

  const Rot3 R(local.expmap());  // nRb: rotation of body in nav frame

  const Vector3 alpha_tangent =  // angular acceleration mapped back to tangent space
      local.applyInvDexp(alpha_body, A ? &alpha_tangent_H_theta : 0, C ? &invH : 0);
  // petrlmat: alpha_tangent_H_theta = 0.5 * skewSymmetric(alpha_body)

  const Vector3 a_nav = R * a_body;
  const double  dt22  = 0.5 * dt * dt;

  Vector12 preintegratedPlus;
  preintegratedPlus <<  // new preintegrated vector:
      theta + omega * dt + alpha_tangent * dt22,  // theta
      position + velocity * dt + a_nav * dt22,        // position
      velocity + a_nav * dt,                          // velocity
      omega + alpha_tangent * dt;                     // omega

  if (A) {
    // Exact derivative of R*a with respect to theta:
    const Matrix3 a_nav_H_theta = R.matrix() * skewSymmetric(-a_body) * local.dexp();
    // Exact derivative of R*alpha with respect to theta:

    A->setIdentity();
    A->block<3, 3>(0, 0).noalias() += 0.5 * alpha_tangent_H_theta * dt22;  // theta
    A->block<3, 3>(3, 0) = a_nav_H_theta * dt22;    // position wrpt theta...
    A->block<3, 3>(3, 6) = I_3x3 * dt;              // .. and velocity
    A->block<3, 3>(6, 0) = a_nav_H_theta * dt;      // velocity wrpt theta
    A->block<3, 3>(9, 0) = alpha_tangent_H_theta * dt;  // omega wrpt theta
    A->block<3, 3>(0, 9) = I_3x3 * dt;              // theta wrpt omega
  }
  // Jacobian wrpt a
  if (B) {
    B->block<3, 3>(0, 0) = Z_3x3;
    B->block<3, 3>(3, 0) = R.matrix() * dt22;
    B->block<3, 3>(6, 0) = R.matrix() * dt;
    B->block<3, 3>(9, 0) = Z_3x3;
  }
  // Jacobian wrpt alpha
  if (C) {
    C->block<3, 3>(0, 0) = invH * dt22;
    C->block<3, 3>(3, 0) = Z_3x3;
    C->block<3, 3>(6, 0) = Z_3x3;
    C->block<3, 3>(9, 0) = invH * dt;
  }

  return preintegratedPlus;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ update() */
void MasPreintegration::update(const Vector3& measuredAcc, const Vector3& measuredAlpha, const double dt, Matrix_12x12* A, Matrix_12x3* B,
                                      Matrix_12x3* C) {

  // Correct for bias in the sensor frame
  Vector3 acc;
  if (boost::static_pointer_cast<MasParams>(params())->getUseBiasAccCorrection()) {
    acc = biasHat_.correctLinAcc(measuredAcc);
  } else {
    acc = measuredAcc;
  }

  Vector3 alpha;
  if (boost::static_pointer_cast<MasParams>(params())->getUseBiasAlphaCorrection()) {
    alpha = biasHat_.correctAngAcc(measuredAlpha);  
  } else {
    alpha = measuredAlpha;
  }

  // Do update
  deltaTij_ += dt;
  preintegrated_ = UpdatePreintegrated(acc, alpha, dt, preintegrated_, A, B, C);

  // new_H_biasAcc = new_H_old * old_H_biasAcc + new_H_acc * acc_H_biasAcc
  // where acc_H_biasAcc = -I_3x3, hence
  // new_H_biasAcc = new_H_old * old_H_biasAcc - new_H_acc
  preintegrated_H_biasAcc_ = (*A) * preintegrated_H_biasAcc_ - (*B);

  // new_H_biasOmega = new_H_old * old_H_biasOmega + new_H_omega * omega_H_biasOmega
  // where omega_H_biasOmega = -I_3x3, hence
  // new_H_biasOmega = new_H_old * old_H_biasOmega - new_H_omega
  preintegrated_H_biasAlpha_ = (*A) * preintegrated_H_biasAlpha_ - (*C);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ biasCorrectedDelta() */
Vector12 MasPreintegration::biasCorrectedDelta(const mas_bias::ConstantBias& bias_i, OptionalJacobian<12, 6> H) const {

  // We correct for a change between bias_i and the biasHat_ used to integrate
  // This is a simple linear correction with obvious derivatives
  const mas_bias::ConstantBias biasIncr                  = bias_i - biasHat_;
  Vector12                    biasCorrected             = preintegrated();
  Matrix_12x3                 D_biasCorrected_biasAcc   = Z_12x3;
  Matrix_12x3                 D_biasCorrected_biasAlpha = Z_12x3;
  if (boost::static_pointer_cast<MasParams>(params())->getUseBiasAccCorrection()) {
    biasCorrected += preintegrated_H_biasAcc_ * biasIncr.linAcc();
    if (H) {
      D_biasCorrected_biasAcc << preintegrated_H_biasAcc_;
    }
  }

  if (boost::static_pointer_cast<MasParams>(params())->getUseBiasAlphaCorrection()) {
    biasCorrected += preintegrated_H_biasAlpha_ * biasIncr.angAcc();
    if (H) {
      (*H) << Z_12x3, preintegrated_H_biasAlpha_;
      D_biasCorrected_biasAlpha << preintegrated_H_biasAlpha_;
    }
  }

  if (H) {
    (*H) << D_biasCorrected_biasAcc, D_biasCorrected_biasAlpha;
  }

  return biasCorrected;
}
/*//}*/

/*//{ derivative sugar */
// sugar for derivative blocks
// NOTE(petrlmat): t is not theta but position here!
#define D_R_R(H) (H)->block<3, 3>(0, 0)
#define D_R_t(H) (H)->block<3, 3>(0, 3)
#define D_R_v(H) (H)->block<3, 3>(0, 6)
#define D_R_w(H) (H)->block<3, 3>(0, 9)
#define D_t_R(H) (H)->block<3, 3>(3, 0)
#define D_t_t(H) (H)->block<3, 3>(3, 3)
#define D_t_v(H) (H)->block<3, 3>(3, 6)
#define D_t_w(H) (H)->block<3, 3>(3, 9)
#define D_v_R(H) (H)->block<3, 3>(6, 0)
#define D_v_t(H) (H)->block<3, 3>(6, 3)
#define D_v_v(H) (H)->block<3, 3>(6, 6)
#define D_v_w(H) (H)->block<3, 3>(6, 9)
#define D_w_R(H) (H)->block<3, 3>(9, 0)
#define D_w_t(H) (H)->block<3, 3>(9, 3)
#define D_w_v(H) (H)->block<3, 3>(9, 6)
#define D_w_w(H) (H)->block<3, 3>(9, 9)
/*//}*/

//------------------------------------------------------------------------------

/*//{ Compose() */
Vector12 MasPreintegration::Compose(const Vector12& zeta01, const Vector12& zeta12, double deltaT12, OptionalJacobian<12, 12> H1,
                                           OptionalJacobian<12, 12> H2) {
  const auto t01 = zeta01.segment<3>(0);
  const auto p01 = zeta01.segment<3>(3);
  const auto v01 = zeta01.segment<3>(6);
  const auto w01 = zeta01.segment<3>(9);

  const auto t12 = zeta12.segment<3>(0);
  const auto p12 = zeta12.segment<3>(3);
  const auto v12 = zeta12.segment<3>(6);
  const auto w12 = zeta12.segment<3>(9);

  Matrix3    R01_H_t01, R12_H_t12;
  const Rot3 R01 = Rot3::Expmap(t01, R01_H_t01);
  const Rot3 R12 = Rot3::Expmap(t12, R12_H_t12);

  // NOTE(petrlmat): on-manifold composition of rotation
  Matrix3    R02_H_R01, R02_H_R12;  // NOTE(frank): R02_H_R12 == Identity
  const Rot3 R02 = R01.compose(R12, R02_H_R01, R02_H_R12);

  Matrix3       t02_H_R02;
  Vector12      zeta02;
  const Matrix3 R = R01.matrix();
  /* zeta02 << Rot3::Logmap(R02, t02_H_R02), // theta */
  zeta02 << Rot3::Logmap(R02, t02_H_R02) + w01 * deltaT12,  // theta NOTE(petrlmat): added integration of w01
      p01 + v01 * deltaT12 + R * p12,                       // position
      v01 + R * v12,                                        // velocity
      w01 + R * w12;                                        // omega

  if (H1) {
    H1->setIdentity();
    D_R_R(H1) = t02_H_R02 * R02_H_R01 * R01_H_t01;
    D_R_w(H1) = I_3x3 * deltaT12;  // NOTE(petrlmat): added partial derivative of R wrpt w
    D_t_R(H1) = R * skewSymmetric(-p12) * R01_H_t01;
    D_t_v(H1) = I_3x3 * deltaT12;
    D_v_R(H1) = R * skewSymmetric(-v12) * R01_H_t01;
    D_w_R(H1) = R * skewSymmetric(-w12);  // NOTE(petrlmat): added partial derivative of w wrpt R
  }

  if (H2) {
    H2->setZero();
    D_R_R(H2) = t02_H_R02 * R02_H_R12 * R12_H_t12;
    D_t_t(H2) = R;
    D_v_v(H2) = R;
    D_w_w(H2) = R;
  }

  return zeta02;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ mergeWith() */
void MasPreintegration::mergeWith(const MasPreintegration& pim12, Matrix_12x12* H1, Matrix_12x12* H2) {
  if (!matchesParamsWith(pim12)) {
    throw std::domain_error("Cannot merge pre-integrated measurements with different params");
  }

  if (params()->body_P_sensor) {
    throw std::domain_error("Cannot merge pre-integrated measurements with sensor pose yet");
  }

  const double t01 = deltaTij();
  const double t12 = pim12.deltaTij();
  deltaTij_        = t01 + t12;

  const Vector12 zeta01 = preintegrated();
  Vector12       zeta12 = pim12.preintegrated();  // will be modified.

  const mas_bias::ConstantBias bias_incr_for_12 = biasHat() - pim12.biasHat();
  zeta12 += pim12.preintegrated_H_biasAlpha_ * bias_incr_for_12.angAcc()  // TODO petrlmat fix the bias increment
            + pim12.preintegrated_H_biasAcc_ * bias_incr_for_12.linAcc();

  preintegrated_ = MasPreintegration::Compose(zeta01, zeta12, t12, H1, H2);

  preintegrated_H_biasAcc_ = (*H1) * preintegrated_H_biasAcc_ + (*H2) * pim12.preintegrated_H_biasAcc_;

  preintegrated_H_biasAlpha_ = (*H1) * preintegrated_H_biasAlpha_ + (*H2) * pim12.preintegrated_H_biasAlpha_;
}
/*//}*/

//------------------------------------------------------------------------------

}  // namespace gtsam
