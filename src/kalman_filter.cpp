#include "kalman_filter.h"
#include "tools.h" /*CalculateJacobian*/
#indlude <math.h> /*atan*/ 

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_; 
  MatrixXd Ft = F_.transpose(); 
  P_ = F_ * P_ * Ft * Q_; 
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  //TODO lecture 21: the second value in the coordinate vector is an angle. this must be normalizes s.t. the angle is between -pi and pi (subtract 2pi until it is in the range 
  float px = x_(0); 
  float py = x_(1); 
  float vx = x_(2); 
  float vy = x_(3); 
  float c1 = sprt(px*px-py*py);   
 
  // check validity of cooridnates 
  if (px == 0 && py == 0){
    cout << "ERROR! Invalid values for px and py: both cannot be zero.";
  }
  // Calculate Jacobian H
  MatrixXd Hj = CalculateJacobian(&x_);
  MatrixXd Hjt = Hj.transpose();
    
  // Convert cartesian coordinates to polar coordinates   
  VectorXd z_pred(3); 
  z_pred << c1, 
      atan(py/px),
      (px*vx+py*vy)/c1; 
  
  y = z - z_pred; 
  S = Hj * P_ * Hjt + R_; 
  Si = S.inverse(); 
  K = P_ * Hjt * Si;
  x_ = x_ + (K*y); 
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
