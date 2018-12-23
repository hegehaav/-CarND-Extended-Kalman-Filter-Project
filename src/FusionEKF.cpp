#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "kalman_filter.h"
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
  
  // Initial state covariance matrix 
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1,0,0,0,
             0,1,0,0,
             0,0,1000,0,
             0,0,0,1000; 
  
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

    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      // get polar coordinates from data 
      double rho =     measurement_pack.raw_measurements_[0];
      double phi =     measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2]; 
      
      // convert polar coordinates to cartesian coordinates 
      double px = rho*cos(phi); 
      if (px < 0.0001) {
          px = 0.0001;      
      }
      
      double py = rho*sin(phi); 
      if (py < 0.0001) {
         py = 0.0001;      
      }
      
      double vx = rho_dot*cos(phi); 
      double vy = rho_dot*sin(phi); 

      // initilize state 
      ekf_.x_ << px, py, vx, vy;
    }
  
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0; 
    }
    // update timestamp
    previous_timestamp_ = measurement_pack.timestamp_; 

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  /**
   * Prediction
   */
  
  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Initial transistion matrix 
  ekf_.F_ = MatrixXd(4,4);
  // update the state transition matrix F according to the new elapsed time.
  ekf_.F_ << 1,0,dt,0,
             0,1,0,dt,
             0,0,1,0,
             0,0,0,1;
  
  // time variables 
  double dt_2 = dt*dt; 
  double dt_3 = dt_2*dt; 
  double dt_4 = dt_3*dt; 
  
  // acceleration noise components 
  double noise_ax = 9; 
  double noise_ay = 9;
  
  // update the the process noise covariance matrix.
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  // predict 
  ekf_.Predict();
  cout << "Predict"<<endl; 

  /**
   * Update:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar Update:
    
    // measurement covariance matrix - radar
    MatrixXd R_radar_ = MatrixXd(3,3); 
    R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

    Tools tools; 
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;    
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    cout<<"Update Radar"<<endl;

  } else {
    // Laser Update
    //measurement covariance matrix - laser
    MatrixXd R_laser_ = MatrixXd(2,2);
    R_laser_ << 0.0225, 0,
              0, 0.0225;
    // measurement matrix - laser 
    MatrixXd H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0; 
    ekf_.H_ = H_laser_; 
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    cout << "Update Laser"<<endl; 
  }
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
