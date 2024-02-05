#ifndef UTILITY_H
#define UTILITY_H

#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/Quaternion.h>

#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/geometry/Quaternion.h>

#include "mas_factor/mas_bias.h"

namespace gtsam_playground
{

/*//{ fromMsg() */
gtsam::Vector3 fromMsg(const geometry_msgs::Vector3& vec_in) {
  gtsam::Vector3 vec_out;
  vec_out[0] = vec_in.x;
  vec_out[1] = vec_in.y;
  vec_out[2] = vec_in.z;
  return vec_out;
}

gtsam::Rot3 fromMsg(const geometry_msgs::Quaternion& q_in) {
  return gtsam::Rot3(gtsam::Quaternion(q_in.w, q_in.x, q_in.y, q_in.z));
}
/*//}*/

/* applyOffset //{ */
gtsam::Pose3 applyOffset(const gtsam::Pose3& msg, const gtsam::Pose3& pose_offset) {

  const gtsam::Point3 pos = pose_offset.rotation().inverse() * (msg.translation() - pose_offset.translation());
  const gtsam::Rot3   rot = pose_offset.rotation().inverse() * msg.rotation();

  return gtsam::Pose3(rot, pos);
}

//}

/*//{ saveOutput() */
void saveOutput(const std::unique_ptr<std::ofstream>& file, const double& t, const gtsam::Pose3& pose, const gtsam::Vector3& lin_vel,
                const gtsam::Vector3& ang_vel) {
  std::stringstream ss;
  ss << t << "," << pose.translation()[0] << "," << pose.translation()[1] << "," << pose.translation()[2] << "," << lin_vel[0] << "," << lin_vel[1] << ","
     << lin_vel[2] << "," << pose.rotation().rpy()[0] << "," << pose.rotation().rpy()[1] << "," << pose.rotation().rpy()[2] << "," << ang_vel[0] << ","
     << ang_vel[1] << "," << ang_vel[2] << endl;
  *file << ss.str();
  file->flush();
}
/*//}*/

/*//{ savePreintegration() */
void savePreintegration(const std::unique_ptr<std::ofstream>& file, const double& deltaT, const gtsam::Vector3& last_acc, const gtsam::Vector3& last_alpha,
                        const gtsam::Vector3& deltaP, const gtsam::Vector3& deltaV, const gtsam::Vector3& deltaR, const gtsam::Vector3& deltaW,
                        const gtsam::mas_bias::ConstantBias& bias) {
  std::stringstream ss;
  ss << deltaT << "," << last_acc[0] << "," << last_acc[1] << "," << last_acc[2] << "," << last_alpha[0] << "," << last_alpha[1] << "," << last_alpha[2] << ","
     << deltaP[0] << "," << deltaP[1] << "," << deltaP[2] << "," << deltaV[0] << "," << deltaV[1] << "," << deltaV[2] << "," << deltaR[0] << "," << deltaR[1]
     << "," << deltaR[2] << "," << deltaW[0] << "," << deltaW[1] << "," << deltaW[2] << "," << bias.linAcc()[0] << "," << bias.linAcc()[1] << ","
     << bias.linAcc()[2] << "," << bias.angAcc()[0] << "," << bias.angAcc()[1] << "," << bias.angAcc()[2] << endl;
  *file << ss.str();
  file->flush();
}
/*//}*/

/*//{ saveAccAlpha() */
void saveAccAlpha(const std::unique_ptr<std::ofstream>& file, const double& deltaT, const gtsam::Vector3& last_acc, const gtsam::Vector3& last_alpha,
                  const gtsam::mas_bias::ConstantBias& bias, const gtsam::Vector& mas, const double last_acc_norm) {
  std::stringstream ss;
  ss << deltaT << "," << last_acc[0] << "," << last_acc[1] << "," << last_acc[2] << "," << last_alpha[0] << "," << last_alpha[1] << "," << last_alpha[2] << ","
     << bias.linAcc()[0] << "," << bias.linAcc()[1] << "," << bias.linAcc()[2] << "," << bias.angAcc()[0] << "," << bias.angAcc()[1] << "," << bias.angAcc()[2]
     << "," << mas[0] << "," << mas[1] << "," << mas[2] << "," << mas[3] << "," << last_acc_norm << endl;
  *file << ss.str();
  file->flush();
}
/*//}*/

/*//{ saveJacobian() */
void saveJacobian(const std::unique_ptr<std::ofstream>& file, const gtsam::Matrix& H, const double iter) {
  std::stringstream ss;
  ss << iter << ":" << H.rows() << ',' << H.cols() << endl << H << endl;
  *file << ss.str();
  file->flush();
}
/*//}*/

}  // namespace gtsam_playground

#endif
