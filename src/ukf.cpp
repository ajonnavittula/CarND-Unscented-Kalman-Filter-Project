#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, first measurement will be recorded

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  // If not initialized, record first measurement

  if (!is_initialized_)
  {
    x_ << 1, 1, 1, 1, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      
      // Laser gives position information only

      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0,
            0,
            0; 
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {

      // Convert from polar co-ordinates to cartesian co-ordinates

      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[3];

      x_ << rho * cos(phi),
            rho * sin(phi),
            rho_dot,
            phi,
            0;
    }

    // Set initial probabilities based on course values

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1000, 0, 0,
          0, 0, 0, 1000, 0,
          0, 0, 0, 0, 1;

    is_initialized_ = true;

    return;
  }




}

// Used implementation from course as a baseline

MatrixXd UKF::GenerateSignmaPoints() {

  // Get state dimension and set lambda

  int n_x = x.size();

  double lambda = 3 - n_x;

  // Get P inverse using Cholesky decomposition

  MatrixXd A = P.llt().matrixL();

  // Create matrix to store sigma points
  
  MatrixXd Xsig = MatrixXd(n_x, 2*n_x+1);

  Xsig.col(0) << x_;

  MatrixXd sigmas_ = MatrixXd(n_x, n_x);
  sigmas_ = sqrt(lambda + n_x)*A;

  for (int i=1; i<n_x+1; i++)
  {
    Xsig.col(i) << x + sigmas_.col(i - 1);

    Xsig.col(n_x + i) << x - sigmas_.col(i + n_x - 1);
  }

  return Xsig;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}