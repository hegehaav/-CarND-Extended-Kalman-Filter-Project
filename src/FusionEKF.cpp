#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h> /*tan*/
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  /* TODOs: 
     - should initial value of Hj be set 
     - should x_ be initialized 
     - should Q be initialized 
*/
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4); 
  P_ = MatrixXd(4,4); 
  F_ = MatrixXd(4,4); 
  Q_ = MatrixXd(4,4); 

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // measurement matrix - laser 
  H_laser << 1, 0, 0, 0,
             0, 1, 0, 0; 
  
  // state covariance matrix 
  P_ << 1,0,0,0,
        0,1,0,0,
        0,0,1000,0,
        0,0,0,1000; 
  
  // initial transistion matrix 
  F_ << 1,0,1,0,
        0,1,0,1,
        0,0,1,0,
        0,0,0,1;
  // Initial process covariance matrix 
  Q_ << 0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0;
          
  /**
   * TODO: Finish initializing the FusionEKF.
   * todo: Set the process and measurement noises
       
   */


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * todo: Initialize the state ekf_.x_ with the first measurement.
     * todo: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    KalmanFilter ekf_;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.Q_ = Q_;
    ekf_.F_ = F_; 
    ekf_.P_ = P_; 

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // todo: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      
      // get polar coordinates from data 
      float rho = mearuement_pack.raw_measurements_[1];
      float phi = mearuement_pack.raw_measurements_[2]; 
      
      // convert polar coordinates to cartesian coordinates 
      float py = (tan(phi)*rho)/sqrt(1 + (tan(phi)*tan(phi)));
      float px = sqrt(rho*rho - py*py);
      
      // initilize state 
      ekf_.R_ = R_radar; 
      ekf_.x_ << px,
                 py,
                 0,
                 0;
      ekf_.H_ = CalculateJacobian(&ekf_.x_);
      // update timestamp
      previous_timestamp = measurement_pack.timestamp_; 
      
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // todo: Initialize state.

      // initialize state 
      // initialize radar covariance matrix 
      ekf_.R_ = R_laser; 
      ekf_.H_ = H_laser; 
      ekf_.x_ << mearuement_pack.raw_measurements_[1],
                 mearuement_pack.raw_measurements_[2],
                 0,
                 0; 
      // update timestamp
      previous_timestamp = measurement_pack.timestamp_; 

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  /**
   * Prediction
   */

  /**
   * todo: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * todo: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // time variables 
  float dt_2 = dt*dt; 
  float dt_3 = dt_2*dt; 
  float dt_4 = dt_3*dt; 
  
  // acceleration noise components 
  noise_ax = 9; 
  noise_ay = 9;
  
  // update the state transition matrix F according to the new elapsed time.
  ekf_.F_(0,2) = dt; 
  ekf_.F_(1,3) = dt; 
  
  // update the the process noise covariance matrix.
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;



  ekf_.Predict();

  /**
   * Update
   */

  /**
   * todo:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // todo: Radar updates
    ekf_.Update(measurement_pack.raw_measurements_);

  } else {
    // todo: Laser updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
