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
  
  // Augmented state vector
  x_aug_ = VectorXd(7);

  // Augmented process covariance matrix
  P_aug_ = MatrixXd::Zero(7, 7);

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

  n_x = x.size();

  int n_aug = 7;

  double lambda = 3 - n_aug;

  x_aug_.head(5) = x;
  x_aug_.tail(2) << 0,
                    0;

  MatrixXd Q = MatrixXd(2,2);

  Q << pow(std_a,2) ,0,
       0, pow(std_yawdd_,2);

  P_aug_.topLeftCorner(n_x,n_x) = P_;
  P_aug_.bottomRightCorner(2,2) = Q;

  // Get sqrt of P using Cholesky decomposition

  MatrixXd A = P_aug_.llt().matrixL();

  // Create matrix to store sigma points
  
  MatrixXd Xsig = MatrixXd(n_aug, 2*n_aug+1);

  Xsig.col(0) << x_aug_;

  MatrixXd sigmas_ = MatrixXd(n_aug, n_aug);
  sigmas_ = sqrt(lambda + n_aug)*A;

  for (int i=1; i<n_aug+1; i++)
  {
    Xsig.col(i) << x_aug_ + sigmas_.col(i - 1);

    Xsig.col(n_aug + i) << x_aug_ - sigmas_.col(i + n_aug - 1);
  }

  return Xsig;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  MatrixXd Xsig_aug = UKF::GenerateSignmaPoints();
  
  Xsig_pred_ = MatrixXd(n_x, 2 * n_aug + 1);

  VectorXd x_k = VectorXd(n_aug);
  VectorXd v = VectorXd(n_x);
  VectorXd a = VectorXd(n_x);
  
  float v_k, psi_k_dot, psi_k, nu_a_k, nu_psi_k; 
  
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
      x_k = Xsig_aug.col(i);
      v_k = x_k(2);
      psi_k = x_k(3);
      psi_k_dot = x_k(4);
      nu_a_k = x_k(5);
      nu_psi_k = x_k(6);
      
      a << 0.5 * nu_a_k * cos(psi_k) * pow(delta_t,2),
           0.5 * nu_a_k * sin(psi_k) * pow(delta_t,2),
           delta_t * nu_a_k,
           0.5 * pow(delta_t,2) * nu_psi_k,
           delta_t * nu_psi_k;
      
      if (psi_k_dot != 0)
      {
          v << (sin(psi_k + psi_k_dot * delta_t) - sin(psi_k)) * v_k / psi_k_dot,
               (-cos(psi_k + psi_k_dot * delta_t) + cos(psi_k)) * v_k / psi_k_dot,
               0,
               psi_k_dot * delta_t,
               0;
      }
      else
      {
          v << v_k * cos(psi_k) * delta_t,
               v_k * sin(psi_k) * delta_t,
               0,
               psi_k_dot * delta_t,
               0;
      }
      
      Xsig_pred_.col(i) = x_k.head(5) + v + a;
  }

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