#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl; 

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Initial state covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /* Initialization */
  if (!is_initialized_) {
    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      // get polar coordinates from data 
      double rho = measurement_pack.raw_measurements_[0]; 
      double phi = measurement_pack.raw_measurements_[1]; 
      double rho_dot = measurement_pack.raw_measurements_[2]; 
      
      // convert polar coordinates to cartesian coordinates 
      double px = rho * cos(phi);
      if ( px < 0.0001 ) {
        px = 0.0001;
      }
      double py = rho * sin(phi);
      if ( py < 0.0001 ) {
        py = 0.0001;
      }
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      ekf_.x_ << px, py, vx , vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // update timestamp
    previous_timestamp_ = measurement_pack.timestamp_ ;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /* Prediction */

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
   double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

  // state transition matrix 
  ekf_.F_ = MatrixXd(4, 4);
  // update the state transition matrix F according to the new elapsed time.
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  double dt_2 = dt * dt; //dt^2
  double dt_3 = dt_2 * dt; //dt^3
  double dt_4 = dt_3 * dt; //dt^4
  
  double noise_ax = 9.0;
  double noise_ay = 9.0;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4/4 * noise_ax, 0, dt_3/2 * noise_ax, 0,
             0, dt_4/4 * noise_ay, 0, dt_3/2 * noise_ay,
             dt_3/2 * noise_ax, 0, dt_2 * noise_ax, 0,
             0, dt_3/2 * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /* Update: update based on sensor type */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else {
    // Laser update
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}