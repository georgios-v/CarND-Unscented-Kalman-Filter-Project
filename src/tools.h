#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Tools {
public:
	/**
	 * Constructor.
	 */
	Tools();

	/**
	 * Destructor.
	 */
	virtual ~Tools();

	/**
	 * A helper method to calculate RMSE.
	 */
	static VectorXd CalculateRMSE(const std::vector<VectorXd> &estimations, const std::vector<VectorXd> &ground_truth);
	
	static double CalculeNIS(const VectorXd &Dz, const MatrixXd &Si);
	
	static void PrintNIS(MeasurementPackage::SensorType sensor, std::vector<double>& values);
	
	static void PrintGraph(const std::vector<double> &data, const std::string fileName, const double reference, 
		const std::string title, const std::string xTitle, const std::string yTitle);
};

#endif /* TOOLS_H_ */