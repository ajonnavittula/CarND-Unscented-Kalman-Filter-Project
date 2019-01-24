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

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

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

    P_ = MatrixXd::Identity(7, 7);

    is_initialized_ = true;

    return;
  }




}

// Used implementation from course as a baseline

MatrixXd UKF::GenerateSigmaPoints() {

  // Get state dimension and set lambda

  x_aug_.head(5) = x_;
  x_aug_.tail(2) << 0,
                    0;

  MatrixXd Q = MatrixXd(2,2);

  Q << pow(std_a_,2) ,0,
       0, pow(std_yawdd_,2);

  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_.bottomRightCorner(2,2) = Q;

  // Get sqrt of P using Cholesky decomposition

  MatrixXd A = P_aug_.llt().matrixL();

  // Create matrix to store sigma points
  
  MatrixXd Xsig = MatrixXd(n_aug_, 2*n_aug_+1);

  Xsig.col(0) << x_aug_;

  MatrixXd sigmas_ = MatrixXd(n_aug_, n_aug_);
  sigmas_ = sqrt(lambda_ + n_aug_)*A;

  for (int i=1; i<n_aug_+1; i++)
  {
    Xsig.col(i) << x_aug_ + sigmas_.col(i - 1);

    Xsig.col(n_aug_ + i) << x_aug_ - sigmas_.col(i + n_aug_ - 1);
  }

  return Xsig;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  MatrixXd Xsig_aug = UKF::GenerateSigmaPoints();
  
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  VectorXd x_k = VectorXd(n_aug_);
  VectorXd v = VectorXd(n_x_);
  VectorXd a = VectorXd(n_x_);
  
  float v_k, psi_k_dot, psi_k, nu_a_k, nu_psi_k; 
  
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
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

      weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // Update state vector with predictions

  x_.fill(0.);

  for(int i = 0; i < n_x; i++)
  {
    x_(i) = weights_.dot(Xsig_pred_.row(i));
  }

  P.fill(0.);

  for(int i = 0; i < 2 * n_aug + 1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
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
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z); 

  float p_x, p_y, v, psi, psi_dot, rho, rho_dot, phi;
  
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
      p_x = Xsig_pred_.col(i)[0];
      p_y = Xsig_pred_.col(i)[1];
      v = Xsig_pred_.col(i)[2];
      psi = Xsig_pred_.col(i)[3];
      psi_dot = Xsig_pred_.col(i)[4];
      
      rho = sqrt(p_x * p_x + p_y * p_y);
      phi = atan2(p_y,p_x);
      rho_dot = 1/rho * (p_x * v * cos(psi) + p_y * v * sin(psi));
      
      Zsig.col(i) << rho,
                     phi,
                     rho_dot;
  }
  
  // calculate mean predicted measurement
  
  z_pred.fill(0.);
  for(int i = 0; i < 2 * n_aug + 1; i++)
  {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(3,3);
  R << std_radr * std_radr, 0, 0,
       0, std_radphi * std_radphi, 0,
       0, 0, std_radrd * std_radrd;
       
  S.fill(0.);       
  for (int i =0; i <2 * n_aug + 1; i++)
  {
      S = S + weights(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
  }
  
  S = S + R;

}