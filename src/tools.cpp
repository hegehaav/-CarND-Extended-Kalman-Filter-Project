#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout; 
using std::endl;
Tools::Tools() {}

Tools::~Tools() {}

// This code is more or less the same as explained on the conferences.
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    unsigned int n = estimations.size(); 
    if(n == 0){
      cout << "ERROR! Estimate vector is empty" << endl;
      return rmse;
    }

    if(ground_truth.size() == 0){
      cout << "ERROR!  Ground-truth vector is empty" << endl;
      return rmse;
    }

    if(n != ground_truth.size()){
      cout << "ERROR! Estimation and ground truth are not the same size" << endl;
      return rmse;
    }

    for(unsigned int i=0; i <n; ++i){
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array()*residual.array();
      rmse += residual;
    }
    // calculate the mean 
    rmse = rmse / n;
    // calculate the square root 
    rmse = rmse.array().sqrt();
    // return RMSE 
    return rmse;
}

// This code is more or less the same as explained on the conferences.
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // check state vaøidity 
  if ( x_state.size() != 4 ) {
    cout << "ERROR! The state has the wrong size" << endl;
    return Hj;
  }
    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    double c1 = px*px+py*py;
    double c2 = sqrt(c1);
    double c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "ERROR! x and y cannot both e zero. Initial H_j is returned" << endl;
        return Hj;
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
         py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}