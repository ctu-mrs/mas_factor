#include "full_state.h"

using namespace std;

namespace gtsam
{

#define TIE(R, t, v, w, x)     \
  const Rot3&      R = (x).R_; \
  const Point3&    t = (x).t_; \
  const Velocity3& v = (x).v_; \
  const Velocity3& w = (x).w_;

//------------------------------------------------------------------------------

/*//{ Create() */
FullState FullState::Create(const Rot3& R, const Point3& t, const Velocity3& v, const Velocity3& w, OptionalJacobian<12, 3> H1, OptionalJacobian<12, 3> H2,
                            OptionalJacobian<12, 3> H3, OptionalJacobian<12, 3> H4) {
  if (H1)
    *H1 << I_3x3, Z_3x3, Z_3x3, Z_3x3;
  if (H2)
    *H2 << Z_3x3, R.transpose(), Z_3x3, Z_3x3;
  if (H3)
    *H3 << Z_3x3, Z_3x3, R.transpose(), Z_3x3;
  if (H4)
    *H4 << Z_3x3, Z_3x3, Z_3x3, R.transpose();
  return FullState(R, t, v, w);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ FromPoseVelocity() */
FullState FullState::FromPoseVelocity(const Pose3& pose, const Vector3& lin_vel, const Vector3& ang_vel, OptionalJacobian<12, 6> H1, OptionalJacobian<12, 3> H2,
                                      OptionalJacobian<12, 3> H3) {
  if (H1)
    *H1 << I_3x3, Z_3x3, Z_3x3, I_3x3, Z_3x3, Z_3x3, Z_3x3, Z_3x3;
  if (H2)
    *H2 << Z_3x3, Z_3x3, pose.rotation().transpose(), Z_3x3;
  if (H3)
    *H3 << Z_3x3, Z_3x3, Z_3x3, pose.rotation().transpose();
  return FullState(pose, lin_vel, ang_vel);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ attitude() */
const Rot3& FullState::attitude(OptionalJacobian<3, 12> H) const {
  if (H)
    *H << I_3x3, Z_3x3, Z_3x3, Z_3x3;
  return R_;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ position() */
const Point3& FullState::position(OptionalJacobian<3, 12> H) const {
  if (H)
    *H << Z_3x3, R(), Z_3x3, Z_3x3;
  return t_;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ linVelocity() */
const Vector3& FullState::linVelocity(OptionalJacobian<3, 12> H) const {
  if (H)
    *H << Z_3x3, Z_3x3, R(), Z_3x3;
  return v_;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ angVelocity() */
const Vector3& FullState::angVelocity(OptionalJacobian<3, 12> H) const {
  if (H)
    *H << Z_3x3, Z_3x3, Z_3x3, R();
  return w_;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ bodyLinVelocity() */
Vector3 FullState::bodyLinVelocity(OptionalJacobian<3, 12> H) const {
  const Rot3&    nRb = R_;
  const Vector3& n_v = v_;
  Matrix3        D_bv_nRb;
  Vector3        b_v = nRb.unrotate(n_v, H ? &D_bv_nRb : 0);
  if (H)
    *H << D_bv_nRb, Z_3x3, I_3x3, Z_3x3;
  return b_v;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ bodyAngVelocity() */
Vector3 FullState::bodyAngVelocity(OptionalJacobian<3, 12> H) const {
  const Rot3&    nRb = R_;
  const Vector3& n_w = w_;
  Matrix3        D_bw_nRb;
  Vector3        b_w = nRb.unrotate(n_w, H ? &D_bw_nRb : 0);
  if (H)
    *H << D_bw_nRb, Z_3x3, Z_3x3, I_3x3;
  return b_w;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ matrix() */
Matrix_10x10 FullState::matrix() const {
  Matrix3      R = this->R();
  Matrix_10x10 T;
  T << R, Z_3x3, Z_3x3, t(), Z_3x3, R, Z_3x3, v(), Z_3x3, Z_3x3, R, w(), Vector10::Zero().transpose(), 1.0;
  return T;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ operator<<() */
ostream& operator<<(ostream& os, const FullState& state) {
  os << "R:" << state.attitude();
  os << "p:" << state.position() << endl;
  os << "v:" << Point3(state.linVelocity()) << endl;
  os << "w:" << Point3(state.angVelocity()) << endl;
  return os;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ print() */
void FullState::print(const string& s) const {
  cout << s << *this << endl;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ equals() */
bool FullState::equals(const FullState& other, double tol) const {
  return R_.equals(other.R_, tol) && traits<Point3>::Equals(t_, other.t_, tol) && equal_with_abs_tol(v_, other.v_, tol) &&
         equal_with_abs_tol(w_, other.w_, tol);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ retract() */
FullState FullState::retract(const Vector12&          xi,  //
                             OptionalJacobian<12, 12> H1, OptionalJacobian<12, 12> H2) const {
  TIE(nRb, n_t, n_v, n_w, *this);
  Matrix3      D_bRc_xi, D_R_nRb, D_t_nRb, D_v_nRb, D_w_nRb;
  const Rot3   bRc = Rot3::Expmap(dR(xi), H2 ? &D_bRc_xi : 0);
  const Rot3   nRc = nRb.compose(bRc, H1 ? &D_R_nRb : 0);
  const Point3 t   = n_t + nRb.rotate(dP(xi), H1 ? &D_t_nRb : 0);
  const Point3 v   = n_v + nRb.rotate(dV(xi), H1 ? &D_v_nRb : 0);
  const Point3 w   = n_w + nRb.rotate(dW(xi), H1 ? &D_w_nRb : 0);
  /* const Point3 w = n_w + dW(xi); */
  if (H1) {
    *H1 << D_R_nRb, Z_3x3, Z_3x3, Z_3x3,  //
        // We pre-multiply with nRc' to account for FullState::Create
        // Then we make use of the identity nRc' * nRb = bRc'
        nRc.transpose() * D_t_nRb, bRc.transpose(), Z_3x3, Z_3x3,
        // Similar reasoning for v:
        nRc.transpose() * D_v_nRb, Z_3x3, bRc.transpose(), Z_3x3, 
        nRc.transpose() * D_w_nRb, Z_3x3, Z_3x3, bRc.transpose();
    /* Z_3x3, Z_3x3, Z_3x3, bRc.transpose(); */
  }
  if (H2) {
    *H2 << D_bRc_xi, Z_3x3, Z_3x3, Z_3x3,      //
        Z_3x3, bRc.transpose(), Z_3x3, Z_3x3,  //
        Z_3x3, Z_3x3, bRc.transpose(), Z_3x3,  //
        Z_3x3, Z_3x3, Z_3x3, bRc.transpose();
  }
  return FullState(nRc, t, v, w);
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ localCoordinates() */
Vector12 FullState::localCoordinates(const FullState&         g,  //
                                     OptionalJacobian<12, 12> H1, OptionalJacobian<12, 12> H2) const {
  Matrix3      D_dR_R, D_dt_R, D_dv_R, D_dw_R;
  const Rot3   dR = R_.between(g.R_, H1 ? &D_dR_R : 0);
  const Point3 dt = R_.unrotate(g.t_ - t_, H1 ? &D_dt_R : 0);
  const Vector dv = R_.unrotate(g.v_ - v_, H1 ? &D_dv_R : 0);
  const Vector dw = R_.unrotate(g.w_ - w_, H1 ? &D_dw_R : 0);

  Vector12 xi;
  Matrix3  D_xi_R;
  xi << Rot3::Logmap(dR, (H1 || H2) ? &D_xi_R : 0), dt, dv, dw;
  if (H1) {
    *H1 << D_xi_R * D_dR_R, Z_3x3, Z_3x3, Z_3x3,  //
        D_dt_R, -I_3x3, Z_3x3, Z_3x3,             //
        D_dv_R, Z_3x3, -I_3x3, Z_3x3,             //
        D_dw_R, Z_3x3, Z_3x3, -I_3x3;
  }
  if (H2) {
    *H2 << D_xi_R, Z_3x3, Z_3x3, Z_3x3,    //
        Z_3x3, dR.matrix(), Z_3x3, Z_3x3,  //
        Z_3x3, Z_3x3, dR.matrix(), Z_3x3,  //
        Z_3x3, Z_3x3, Z_3x3, dR.matrix();
  }
  return xi;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ derivative sugar */
// sugar for derivative blocks
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

/*//{ update() */
FullState FullState::update(const Vector3& b_acceleration, const Vector3& b_alpha, const double dt, OptionalJacobian<12, 12> F, OptionalJacobian<12, 3> G1,
                            OptionalJacobian<12, 3> G2) const {

  Vector12    xi;
  Matrix_3x12 D_xiP_state, D_xiR_state;
  Vector3     b_v  = bodyLinVelocity(F ? &D_xiP_state : 0);
  Vector3     b_w  = bodyAngVelocity(F ? &D_xiR_state : 0);
  double      dt22 = 0.5 * dt * dt;

  // Integrate on tangent space.
  /* dR(xi) << dt * b_w; */
  dR(xi) << dt * b_w + dt22 * b_alpha;  // NOTE(petrlmat): both omega and alpha on the same tangent space
  dP(xi) << dt * b_v + dt22 * b_acceleration;
  dV(xi) << dt * b_acceleration;
  dW(xi) << dt * b_alpha;

  // Bring back to manifold
  Matrix_12x12 D_newState_xi;
  FullState    newState = retract(xi, F, G1 || G2 ? &D_newState_xi : 0);

  // Derivative wrt state is computed by retract directly
  // However, as dP(xi) also depends on state, we need to add that contribution
  if (F) {
    F->middleRows<3>(3) += dt * D_t_t(F) * D_xiP_state;
    F->topRows<3>() += dt * D_R_R(F) * D_xiR_state;  // NOTE(petrlmat): dR(xi) also depends on state (b_w)
  }
  // derivative wrt acceleration
  if (G1) {
    // D_newState_dPxi = D_newState_xi.middleCols<3>(3)
    // D_dPxi_acc = dt22 * I_3x3
    // D_newState_dVxi = D_newState_xi.rightCols<3>()
    // D_dVxi_acc = dt * I_3x3
    // *G2 = D_newState_acc = D_newState_dPxi * D_dPxi_acc + D_newState_dVxi * D_dVxi_acc
    *G1 = D_newState_xi.middleCols<3>(3) * dt22 + D_newState_xi.middleCols<3>(6) * dt;
  }
  // derivative wrt alpha
  if (G2) {
    // D_newState_dRxi = D_newState_xi.leftCols<3>()
    // D_dRxi_omega = dt * I_3x3
    // *G1 = D_newState_omega = D_newState_dRxi * D_dRxi_omega
    *G2 = D_newState_xi.leftCols<3>() * dt22
          + D_newState_xi.rightCols<3>() * dt;
  }
  return newState;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ coriolis() */
// TODO(petrlmat) how to handle?
Vector12 FullState::coriolis(double dt, const Vector3& omega, bool secondOrder, OptionalJacobian<12, 12> H) const {
  TIE(nRb, n_t, n_v, n_w, *this);
  cout << "calculating coriolis (not yet handled!)" << endl;
  const double  dt2             = dt * dt;
  const Vector3 omega_cross_vel = omega.cross(n_v);

  Vector12 xi;
  Matrix3  D_dP_R;
  dR(xi) << nRb.unrotate((-dt) * omega, H ? &D_dP_R : 0);
  dP(xi) << ((-dt2) * omega_cross_vel);  // NOTE(luca): we got rid of the 2 wrt INS paper
  dV(xi) << ((-2.0 * dt) * omega_cross_vel);
  if (secondOrder) {
    const Vector3 omega_cross2_t = omega.cross(omega.cross(n_t));
    dP(xi) -= (0.5 * dt2) * omega_cross2_t;
    dV(xi) -= dt * omega_cross2_t;
  }
  if (H) {
    H->setZero();
    const Matrix3 Omega         = skewSymmetric(omega);
    const Matrix3 D_cross_state = Omega * R();
    H->setZero();
    D_R_R(H) << D_dP_R;
    D_t_v(H) << (-dt2) * D_cross_state;
    D_v_v(H) << (-2.0 * dt) * D_cross_state;
    if (secondOrder) {
      const Matrix3 D_cross2_state = Omega * D_cross_state;
      D_t_t(H) -= (0.5 * dt2) * D_cross2_state;
      D_v_t(H) -= dt * D_cross2_state;
    }
  }
  return xi;
}
/*//}*/

//------------------------------------------------------------------------------

/*//{ correctPIM() */
Vector12 FullState::correctPIM(const Vector12& pim, double dt, const Vector3& n_gravity, const boost::optional<Vector3>& omegaCoriolis,
                               bool use2ndOrderCoriolis, OptionalJacobian<12, 12> H1, OptionalJacobian<12, 12> H2) const {
  const Rot3&      nRb  = R_;
  const Velocity3& n_v  = v_;  // derivative is Ri !
  const Velocity3& n_w  = w_;
  const double     dt22 = 0.5 * dt * dt;

  Vector12 xi;
  Matrix3  D_dR_Ri1, D_dP_Ri1, D_dP_Ri2, D_dP_nv, D_dR_nw, D_dV_Ri;
  // NOTE(petrlmat): unrotate: H1 derivative wrt R, H2 wrt p
  /* dR(xi) = dR(pim); */
  dR(xi) = dR(pim) + dt * nRb.unrotate(n_w, H1 ? &D_dR_Ri1 : 0, H2 ? &D_dR_nw : 0);
  /* + dt * nRb.unrotate(n_w, H1 ? &D_dR_Ri1 : 0, &D_dR_nw); */

  dP(xi) = dP(pim)
           /* + dt * nRb.unrotate(n_v, H1 ? &D_dP_Ri1 : 0, &D_dP_nv) */
           + dt * nRb.unrotate(n_v, H1 ? &D_dP_Ri1 : 0, H2 ? &D_dP_nv : 0) + dt22 * nRb.unrotate(n_gravity, H1 ? &D_dP_Ri2 : 0);
           /* + dt * nRb.unrotate(n_v, H1 ? &D_dP_Ri1 : 0, H2 ? &D_dP_nv : 0); */
  dV(xi) = dV(pim) + dt * nRb.unrotate(n_gravity, H1 ? &D_dV_Ri : 0);
  /* dV(xi) = dV(pim); */
  dW(xi) = dW(pim);  // NOTE(petrlmat): no correction to omega so far

  if (omegaCoriolis) {
    // NOTE(petrlmat) not yet handled
    if (H1) {
      H1->setZero();  // NOTE(petrlmat) we are not handling the coriolis yet, so we must at least initialize H1
    }

    /* xi += coriolis(dt, *omegaCoriolis, use2ndOrderCoriolis, H1); */
  }

  if (H1 || H2) {
    Matrix3 Ri = nRb.matrix();

    if (H1) {  // NOTE(petrlmat): should be correct now
      if (!omegaCoriolis)
        H1->setZero();  // if coriolis H1 is already initialized
      D_R_R(H1) += dt * D_dR_Ri1;
      D_R_w(H1) += dt * D_dR_nw * Ri;
      D_t_R(H1) += dt * D_dP_Ri1 + dt22 * D_dP_Ri2;
      /* D_t_R(H1) += dt * D_dP_Ri1; */
      D_t_v(H1) += dt * D_dP_nv * Ri;
      D_v_R(H1) += dt * D_dV_Ri;
      /* cout << "correctPIM H: " << endl << *H1 << endl; */
      /* if (omegaCoriolis) { */
      /* cout << "omegaCoriolis: " << *omegaCoriolis << endl; */
      /* } */
      /* cout << "dt: " << dt << endl; */
      /* cout << "Ri: " << endl << Ri << endl; */
      /* cout << "D_dR_Ri1: " << endl << D_dR_Ri1 << endl; */
      /* cout << "D_dR_nw: " << endl << D_dR_nw << endl; */
      /* cout << "D_dP_Ri1: " << endl << D_dP_Ri1 << endl; */
      /* cout << "D_dP_Ri2: " << endl << D_dP_Ri2 << endl; */
      /* cout << "D_dP_nv: " << endl << D_dP_nv << endl; */
      /* cout << "D_dV_Ri: " << endl << D_dV_Ri << endl; */
    }

    if (H2) {
      H2->setIdentity();  // NOTE(petrlmat) we are not handling the coriolis yet, so we must at least initialize H2
    }
  }

  return xi;
}
/*//}*/

//------------------------------------------------------------------------------

}  // namespace gtsam
