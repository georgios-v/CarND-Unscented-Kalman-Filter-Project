#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>
#include <vector>

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
	std_a_ = 0.6;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = M_PI/6;

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

	previous_timestamp_ = 0;

	initialized_ = false;

	///* State dimension
	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	///* Lidar measurement state dimension
	n_z_l_ = 2;

	///* Radar measurement state dimension
	n_z_r_ = 3;

	///* Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred_.fill(0.0);

	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug_.fill(0.0);

	//set weights
	weights_ = VectorXd(2 * n_aug_ + 1);
	for (int i = 0; i < 2 * n_aug_ + 1; ++i)
		weights_[i] = ((i != 0) * 0.5 + (i == 0) * lambda_) / (lambda_ + n_aug_);

	///* Store NIS values
	nis_lidar_ = vector<double>();
	nis_radar_ = vector<double>();

	///* Use Gnuplot to print the NIS graph
	print_graph_ = false;
}

UKF::~UKF() {
}

void UKF::GenerateAugmentedSigmaPoints() {
	lambda_ = 3 - n_aug_;
	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create augmented mean state
	x_aug << x_, 0, 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug.bottomRightCorner(2, 2) << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();

	//create augmented sigma points
	MatrixXd tmp = sqrt(lambda_ + n_aug_) * A;

	Xsig_aug_.col(0) = x_aug;
	for (int i = 0; i < n_aug_; ++i) {
		Xsig_aug_.col(i + 1) = x_aug + tmp.col(i);
		Xsig_aug_.col(i + 1 + n_aug_) = x_aug - tmp.col(i);
	}

	//print result
	std::cout << "Xsig_aug = " << std::endl << Xsig_aug_ << std::endl;
}

void UKF::SigmaPointPrediction(double delta_t) {
	double delta_t2 = delta_t * delta_t / 2.0;

	VectorXd P_pred(n_x_), Q_pred(n_x_);
	//predict sigma points
	//avoid division by zero
	//write predicted sigma points into right column
	for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
		VectorXd In = Xsig_aug_.col(i);
		if (In[4] == 0) {
			P_pred << In[2] * cos(In[3]) * delta_t,
					In[2] * sin(In[3]) * delta_t,
					0, In[4] * delta_t, 0;
		} else {
			double a = In[2] / In[4];
			double x3 = In[4] * delta_t;
			P_pred << a * (sin(In[3] + x3) - sin(In[3])),
					a * (-cos(In[3] + x3) + cos(In[3])),
					0, x3, 0;
		}
		Q_pred << delta_t2 * cos(In[3]) * In[5],
				delta_t2 * sin(In[3]) * In[5],
				delta_t * In[5],
				delta_t2 * In[6],
				delta_t * In[6];
		Xsig_pred_.col(i) = In.head(n_x_) + P_pred + Q_pred;
	}

	//print result
	std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
}

void UKF::PredictMeanAndCovariance() {
	int n_cols = 2 * n_aug_ + 1;

	//predict state mean
	x_.fill(0.0);
	for (int i = 0; i < n_cols; ++i)
		x_ = x_ + weights_[i] * Xsig_pred_.col(i);

	//predict state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < n_cols; ++i) {
		VectorXd Dx = Xsig_pred_.col(i) - x_;
		while (Dx(3) > M_PI) Dx(3) -= 2. * M_PI;
		while (Dx(3)<-M_PI) Dx(3) += 2. * M_PI;
		P_ = P_ + weights_[i] * Dx * Dx.transpose();
	}

	//print result
	std::cout << "Predicted state" << std::endl;
	std::cout << x_ << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P_ << std::endl;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	if (!initialized_) {
		x_ = VectorXd(n_x_);
		x_ << 1, 1, 1, 1, 0.1;

		// init covariance matrix
		P_ << 0.15, 0, 0, 0, 0,
				0, 0.15, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			/**
			Convert radar from polar to Cartesian coordinates and initialize state.
			 */
			float ro = meas_package.raw_measurements_(0);
			float phi = meas_package.raw_measurements_(1);
			float ro_dot = meas_package.raw_measurements_(2);
			x_ << ro * cos(phi),
					ro * sin(phi),
					ro_dot, 0, 0;
			initialized_ = true;
		} else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			/**
			Initialize state.
			 */
			x_ << meas_package.raw_measurements_, 0, 0, 0;
			initialized_ = true;
		}

		previous_timestamp_ = meas_package.timestamp_;
		return;
	}

	double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = meas_package.timestamp_;

	Prediction(dt);

	if ((meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) ||
			(meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_))
		return;

	Update(meas_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	GenerateAugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();
}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* Zsig_out, MatrixXd* S_out) {
	int n_cols = 2 * n_aug_ + 1;
	//define spreading parameter
	lambda_ = 3 - n_aug_;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_l_, n_cols);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_l_);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_l_, n_z_l_);

	//transform sigma points into measurement space
	Zsig.row(0) = Xsig_pred_.row(0);
	Zsig.row(1) = Xsig_pred_.row(1);

	//calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i = 0; i < n_cols; i++)
		z_pred = z_pred + weights_[i] * Zsig.col(i);

	//calculate innovation covariance matrix S
	S.fill(0.0);
	for (int i = 0; i < n_cols; i++) {
		//residual
		VectorXd Dz = Zsig.col(i) - z_pred;
		//angle normalization
		while (Dz(1) > M_PI) Dz(1) -= 2. * M_PI;
		while (Dz(1)<-M_PI) Dz(1) += 2. * M_PI;
		S = S + weights_[i] * Dz * Dz.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_l_, n_z_l_);
	R << std_laspx_*std_laspx_, 0,
0, std_laspy_*std_laspy_;
	S = S + R;

	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;

	//write result
	*z_out = z_pred;
	*Zsig_out = Zsig;
	*S_out = S;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* Zsig_out, MatrixXd* S_out) {
	int n_cols = 2 * n_aug_ + 1;
	//define spreading parameter
	lambda_ = 3 - n_aug_;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_r_, n_cols);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_r_);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_r_, n_z_r_);

	//transform sigma points into measurement space
	for (int i = 0; i < n_cols; i++) {
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw) * v;
		double v2 = sin(yaw) * v;

		Zsig.col(i) << sqrt(p_x * p_x + p_y * p_y), //r
				atan2(p_y, p_x), //phi
				(p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); //r_dot
	}

	//calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i = 0; i < n_cols; i++)
		z_pred = z_pred + weights_[i] * Zsig.col(i);

	//calculate innovation covariance matrix S
	S.fill(0.0);
	for (int i = 0; i < n_cols; i++) {
		//residual
		VectorXd Dz = Zsig.col(i) - z_pred;

		//angle normalization
		while (Dz(1) > M_PI) Dz(1) -= 2. * M_PI;
		while (Dz(1)<-M_PI) Dz(1) += 2. * M_PI;

		S = S + weights_[i] * Dz * Dz.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_r_, n_z_r_);
	R << std_radr_ * std_radr_, 0, 0,
			0, std_radphi_ * std_radphi_, 0,
			0, 0, std_radrd_ * std_radrd_;
	S = S + R;

	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;

	//write result
	*z_out = z_pred;
	*Zsig_out = Zsig;
	*S_out = S;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(MeasurementPackage meas_package) {
	VectorXd z_pred;
	MatrixXd Zsig;
	MatrixXd S;

	int n_z_ = 0;

	if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		n_z_ = n_z_l_;
		PredictLidarMeasurement(&z_pred, &Zsig, &S);
	} else {
		n_z_ = n_z_r_;
		PredictRadarMeasurement(&z_pred, &Zsig, &S);
	}

	//define spreading parameter
	lambda_ = 3 - n_aug_;

	int n_cols = 2 * n_aug_ + 1;

	//create example vector for incoming radar measurement
	VectorXd& z = meas_package.raw_measurements_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	Tc.fill(0.0);
	for (int i = 0; i < n_cols; ++i) {
		VectorXd Dx = Xsig_pred_.col(i) - x_;
		while (Dx(1) > M_PI) Dx(1) -= 2. * M_PI;
		while (Dx(1)<-M_PI) Dx(1) += 2. * M_PI;

		VectorXd Dz = Zsig.col(i) - z_pred;
		while (Dz(1) > M_PI) Dz(1) -= 2. * M_PI;
		while (Dz(1)<-M_PI) Dz(1) += 2. * M_PI;

		Tc = Tc + weights_[i] * Dx * Dz.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd Si = S.inverse();
	MatrixXd K = Tc * Si;

	//residual
	VectorXd Dz = z - z_pred;
	while (Dz(1) > M_PI) Dz(1) -= 2. * M_PI;
	while (Dz(1)<-M_PI) Dz(1) += 2. * M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * Dz;
	P_ = P_ - K * S * K.transpose();

	//print result
	std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

	//compute NIS
	double nis = Tools::CalculeNIS(Dz, Si);
	std::vector<double>* nis_values;
	if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		nis_values = &nis_lidar_;
		cout << "NIS Lidar: " << nis << endl;
	} else {
		nis_values = &nis_radar_;
		cout << "NIS Radar: " << nis << endl;
	}
	nis_values->push_back(nis);
	//plot Graph
	if (nis_values->size() == 249)
		Tools::PrintNIS(meas_package.sensor_type_, *nis_values);
}
