#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.30;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  time_us_ = 0;
  // State dimensions
  n_x_ = 5;
  // Augmented state dimensions
  n_aug_ = 7;
  // Sigma points lambda
  lambda_ = 3 - n_aug_;
  // initit matrices
  Xsig_pred_ = MatrixXd(n_aug_, n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
      time_us_ = meas_package.timestamp_;
       // initialize state components
      float px; // x pos
      float py; // y pos
      float v = 0.0;
      float theta;
      float theta_ = 0.0;
       // LASER
      if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
           px = meas_package.raw_measurements_[0];
           py = meas_package.raw_measurements_[1];
           theta = atan2(py, px);
           // RADAR
       }
       else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
           float rho = meas_package.raw_measurements_[0];  // range
           float phi = meas_package.raw_measurements_[1];  // bearing
           // Convert from polar to cartesian
           px = rho * cos(phi);
           py = rho * sin(phi);
           theta = atan2(py, px);
       }
      x_ << px, py, v, theta, theta_;
       // Fill with zeros and populate the diagonal
      P_.fill(0.0);
      P_(0, 0) = 1.0;
      P_(1, 1) = 1.0;
      P_(2, 2) = 1.0;
      P_(3, 3) = 1.0;
      P_(4, 4) = 1.0;

       // Done initializing, no need to predict or update
      is_initialized_ = true;
      return;
   }
   /*****************************************************************************
    *  Predict
    ****************************************************************************/
   // time elapsed between the current and previous measurements
   float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
   time_us_ = meas_package.timestamp_;

   Prediction(dt);

   /*****************************************************************************
    *  Update
    ****************************************************************************/
   // LASER update
   if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }
   // RADAR update
   else {
      UpdateRadar(meas_package);
    }
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //Augmemted state vector
  VectorXd x_augmented = VectorXd(n_aug_);
  x_augmented.head(5) = x_;
  x_augmented(5) = 0.0;
  x_augmented(6) = 0.0;
  //Augmented Covariance Matrix
  MatrixXd P_augmented = MatrixXd(n_aug_, n_aug_);
  P_augmented.fill(0.0);
  P_augmented.topLeftCorner(5, 5) = P_;
  P_augmented(5, 5) = std_a_ * std_a_;
  P_augmented(6, 6) = std_yawdd_ * std_yawdd_;
  //Generate Sigma Points Matrix
  MatrixXd AUX = P_augmented.llt().matrixL();
  MatrixXd Xsig_augmented = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_augmented.col(0) = x_augmented;
  for (int i = 0; i < n_aug_; i++) {
      Xsig_augmented.col(i + 1) = x_augmented + sqrt(lambda_ + n_aug_) * AUX.col(i);
      Xsig_augmented.col(i + 1 + n_aug_) =
          x_augmented - sqrt(lambda_ + n_aug_) * AUX.col(i);
  }
  //Predict sigma points
  Xsig_pred_ = pred_sigma_p(Xsig_augmented, delta_t);
  compute_weights_diffs(&x_, &P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;
    MatrixXd H = MatrixXd(2, 5);
    MatrixXd R = MatrixXd(2, 2);
    R << std_laspx_ * std_laspx_, 0,
         0, std_laspy_ * std_laspy_;

    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

    VectorXd z_pred = H * x_;
    VectorXd y = z - z_pred;
    MatrixXd PHt = P_ * H.transpose();
    MatrixXd S = H * PHt + R;
    MatrixXd K = PHt * S.inverse();

    // new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;


    MatrixXd nis = y.transpose() * S.inverse() * y;
    nis_ = nis(0, 0);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    MatrixXd Zsig = meas_sigma_p();
    VectorXd z_pred = pred_meas(Zsig);
    MatrixXd S = meas_cov(Zsig, z_pred);
    VectorXd z = meas_package.raw_measurements_;
    posterior_update(z, z_pred, S, Zsig);

    //NIS
    VectorXd y = z - z_pred;
    MatrixXd nis = y.transpose() * S.inverse() * y;
    nis_ = nis(0, 0);
}


// Theese functions make the calculations with the components for the process
// and keep the code clean and comprehensible;
// inspired and orignally taken from:
// https://github.com/penny4860/CarND-Unscent-Kalman-Filter/blob/master/src/ukf.cpp

// Calculate new Sigma Points
MatrixXd UKF::pred_sigma_p(MatrixXd Xsig_augmented, double delta_t) {
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd x_sigma = Xsig_augmented.col(i);
        VectorXd x_pred(5);

        double px = x_sigma(0);
        double py = x_sigma(1);
        double v = x_sigma(2);
        double theta = x_sigma(3);
        double theta_d = x_sigma(4);
        double n_rad_acc = x_sigma(5);
        double n_theta_acc = x_sigma(6);

        // Avoid division by zero
        if (fabs(theta_d) < 0.001) {
            x_pred(0) = px + v * cos(theta) * delta_t;
            x_pred(1) = py + v * sin(theta) * delta_t;
        }
        else {
            x_pred(0) = px + v / theta_d *
                        (sin(theta + theta_d * delta_t) - sin(theta));
            x_pred(1) = py + v / theta_d *
                        (-cos(theta + theta_d * delta_t) + cos(theta));
        }
        x_pred(2) = v;
        x_pred(3) = theta + theta_d * delta_t;
        x_pred(4) = theta_d;
        // noise
        x_pred(0) += 0.5 * delta_t * delta_t * cos(theta) * n_rad_acc;
        x_pred(1) += 0.5 * delta_t * delta_t * sin(theta) * n_rad_acc;
        x_pred(2) += delta_t * n_rad_acc;
        x_pred(3) += 0.5 * delta_t * delta_t * n_theta_acc;
        x_pred(4) += delta_t * n_theta_acc;

        Xsig_pred.col(i) = x_pred;
    }
    return Xsig_pred;
}
// Operate and predict
void UKF::compute_weights_diffs(VectorXd *x_pred, MatrixXd *P_pred) {

    VectorXd x = VectorXd(n_x_);
    MatrixXd P = MatrixXd(n_x_, n_x_);
    VectorXd weights = VectorXd(2 * n_aug_ + 1);

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        if (i == 0)
            weights(i) = lambda_ / (lambda_ + n_aug_);
        else
            weights(i) = 1 / (2 * (lambda_ + n_aug_));
    }

    // Predict state mean
    x.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        x += weights(i) * Xsig_pred_.col(i);
    }
    // Predict state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        // State differences
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        // Normalize angles
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        P = P + weights(i) * x_diff * x_diff.transpose();
    }
    //Results
    *x_pred = x;
    *P_pred = P;
}

// Sigma Points Matrix in the measurement space
MatrixXd UKF::meas_sigma_p(void) {
    int n_z = 3;
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z, n_z);
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points

        // extract values
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);
        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);  // r
        Zsig(1, i) = atan2(p_y, p_x);              // phi
        Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);  // r_d
    }
    return Zsig;
}

// Calculate sigma weights
VectorXd UKF::sigma_weights(void) {
    VectorXd weights = VectorXd(2 * n_aug_ + 1);
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {
        double weight = 0.5 / (n_aug_ + lambda_);
        weights(i) = weight;
    }
    return weights;
}

//
VectorXd UKF::pred_meas(MatrixXd Zsig) {
    VectorXd weights = sigma_weights();
    // mean predicted measurement
    VectorXd z_pred = VectorXd(3);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    return z_pred;
}

MatrixXd UKF::meas_cov(MatrixXd Zsig, VectorXd z_pred) {
    VectorXd weights = sigma_weights();

    // measurement covariance matrix
    int n_z = 3;
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        S = S + weights(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ * std_radphi_, 0, 0, 0,
        std_radrd_ * std_radrd_;
    S = S + R;
    return S;
}

void UKF::posterior_update(VectorXd z, VectorXd z_pred, MatrixXd S,
                                MatrixXd Zsig) {
    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, 3);
    VectorXd weights = sigma_weights();

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    // residual
    VectorXd z_diff = z - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    // update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
}
