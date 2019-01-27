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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  
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

  // Number of states
  n_x_ = 5;

  // Number of augmented states
  n_aug_ = 7;

  // Predicted sigma point matrix size
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights vector initialization
  weights_ = VectorXd(2 * n_aug_ + 1);

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
      double rho_dot = meas_package.raw_measurements_[2];

      x_ << rho * cos(phi),
            rho * sin(phi),
            rho_dot,
            0,
            0;
    }

    // Set initial probabilities based on course suggestions

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    std::cout << "Initiaization complete" << std::endl;

    return;
  }

  // Calculate time elapsed in seconds
  
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

  // Run prediction step to update state and covariance

  std::cout << "Running prediction" << std::endl;
  Prediction(delta_t);

  // Set current time stamp

  time_us_ = meas_package.timestamp_;
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
   
   std::cout << "Radar measurement received" << std::endl;
   UpdateRadar(meas_package); 
  }
  
  else if ( meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ ) { 
    std::cout << "Lidar complete" << std::endl;
    UpdateLidar(meas_package); 
  }

}


void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  /* Steps required:
     1. Augment state vector with nu
     2. Calculate sigma points
     3. Predict Sigma points using state function
     4. Calculate mean and variance
  */
  /*******************************************************
   **Step 1: Augment state vector and process covariance**
   ******************************************************/

  // Augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);

  // Augmented process covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  // Spreading factor for sigma point calculation
  int lambda = 3 - n_aug_;
  
  // Process noise matrix for augmentation
  MatrixXd Q = MatrixXd(2,2);
  Q << pow(std_a_,2), 0,
       0, pow(std_yawdd_,2);

  // Set x_aug
  x_aug.head(5) = x_;
  x_aug.tail(2) << 0,
                   0;

  // Set P_aug
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q;

  std::cout << "Prediction step 1 complete" << std::endl;
  /**********************************
   **Step 2: Calculate Sigma points**
   *********************************/


  // Sqrt of P using cholesky
  MatrixXd A = P_aug.llt().matrixL();

  MatrixXd sigmas = MatrixXd(n_aug_, n_aug_);
  sigmas = sqrt(lambda + n_aug_) * A;

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  Xsig_aug.col(0) = x_aug;

  for (int i = 1; i < n_aug_ + 1; i++) {
    
    Xsig_aug.col(i) = x_aug + sigmas.col(i - 1);
    Xsig_aug.col(n_aug_ + i) = x_aug - sigmas.col(i - 1);
  }
  
  std::cout << "Prediction step 2 complete" << std::endl;
  /**********************************
   **Step 3: Sigma point prediction**
   *********************************/

  // Some helpful variables
  VectorXd x_k = VectorXd(n_aug_);
  VectorXd v = VectorXd(n_x_);
  VectorXd a = VectorXd(n_x_);

  float v_k, psi_k_dot, psi_k, nu_a_k, nu_psi_k;
  float delta_t_2 = delta_t * delta_t;

  // Since we are already looping, calculate weights for mean calculation
  weights_(0) = lambda / (lambda + n_aug_);


  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    x_k = Xsig_aug.col(i);
    v_k = x_k(2);
    psi_k = x_k(3);
    psi_k_dot = x_k(4);
    nu_a_k = x_k(5);
    nu_psi_k = x_k(6);

    a << 0.5 * cos(psi_k) * nu_a_k * delta_t_2,
         0.5 * sin(psi_k) * nu_a_k * delta_t_2,
         delta_t * nu_a_k,
         0.5 * delta_t_2 * nu_psi_k,
         delta_t * nu_psi_k;

    if (psi_k_dot > 0.01) {
  
      v << v_k / psi_k_dot * (sin(psi_k + psi_k_dot*delta_t) - sin(psi_k)),
      v_k / psi_k_dot * (-cos(psi_k + psi_k_dot*delta_t) + cos(psi_k)),
      0,
      psi_k_dot * delta_t,
      0;
    }

    else {

      v << v_k * cos(psi_k) * delta_t,
           v_k * sin(psi_k) * delta_t,
           0,
           psi_k_dot * delta_t,
           0;
    }

    Xsig_pred_.col(i) = x_k.head(5) + v + a;

    // weight calculation
    if(i > 0) {
      weights_(i) = 0.5 / (lambda + n_aug_);
    }
  }

  std::cout << "Prediction step 3 complete" << std::endl;
  /*************************************************
   **Step 4: mean state and covariance calculation**
   ************************************************/

  x_.fill(0.);
  for (int i = 0; i < n_x_; i++) {

    x_(i) = weights_.dot(Xsig_pred_.row(i));
  }

  P_.fill(0.);
  for (int i = 0; i < 2*n_aug_+1; i++) {

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3) > M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose(); 

    
  }
  std::cout << "Prediction step 4 complete" << std::endl;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  /*
   Using linear equations for update since Lidar is linear 
   Steps required:
   1. Create Laser noise and output matrices
   2. Calculate Kalman gain
   3. Update state
  */

  /**************************************************
   **Step 1: Create laser noise and output matrices**
   *************************************************/
  
  n_z_ = 2;

  //measurement covariance matrix - laser
  MatrixXd R = MatrixXd(n_z_, n_z_);
  MatrixXd H = MatrixXd(n_z_, n_x_);
  VectorXd z = VectorXd(n_z_);

  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1];

  std::cout << "UpdateLidar step 1 complete" << std::endl;
  /*********************************
   **Step 2: Calculate Kalman gain**
   ********************************/

  VectorXd z_pred = H * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  std::cout << "Update Lidar step 2 complete" << std::endl;
  /***************************************
   **Step 3: Update state and covaraince**
   **************************************/
  x_ = x_ + (K * y);

  MatrixXd I =  MatrixXd::Identity(n_x_,n_x_);
  P_ = (I - K * H) * P_;

  std::cout << "Update Lidar step 3 complete" << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: 
   * You can also calculate the radar NIS, if desired.
   */

  /*
   Steps required:
   1. Propagate sigma points through measurement function 
   2. Calculate mean and covariance of measurement sigma points
   3. Update state and covariance using Kalman equations
  */

  /****************************************************************
   **Step 1: Propagate sigma points through measurement function***
   **Step 2: Calcualte mean and covariance using new sigma points**
   **Step 3, Part 1: Calculate variables for Kalman gain***********
   ***************************************************************/
  // Set number of measurement states for Radar
  n_z_ = 3;

  // Incoming measurement vector z
  VectorXd z = VectorXd(3);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1],
       meas_package.raw_measurements_[2];

  // Some helpful variables
  float p_x, p_y, v, psi, rho, rho_dot, phi;

  VectorXd z_pred = VectorXd::Zero(n_z_);

  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  MatrixXd S = MatrixXd(n_z_, n_z_);

  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

  MatrixXd R = MatrixXd(3, 3);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    
    p_x = Xsig_pred_.col(i)[0];
    p_y = Xsig_pred_.col(i)[1];
    v = Xsig_pred_.col(i)[2];
    psi = Xsig_pred_.col(i)[3];

    rho = sqrt( p_x * p_x + p_y * p_y);
    phi = atan2(p_y, p_x);
    rho_dot = 1 / rho * (p_x * v * cos(psi) + p_y * v * sin(psi));

    Zsig.col(i) << rho,
                   phi,
                   rho_dot;

    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;

  S.fill(0.);
  for ( int i = 0; i < 2 * n_aug_ + 1; i++) {

    S = S + weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (x_diff(3) > M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    while (z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  S += R;
  std::cout << "UpdateRadar steps 1,2 and 3 part 1 complete" << std::endl;
  /*****************************************************************
   **Step 3, part 2: Update state and covariance using Kalman gain**
   ****************************************************************/

  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;
  while (z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();

  std::cout << "UpdateRadar step 3 part 2 complete" << std::endl;
}