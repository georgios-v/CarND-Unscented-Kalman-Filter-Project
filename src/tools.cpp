#include <iostream>
#include <vector>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

Tools::Tools() {
}

Tools::~Tools() {
}

VectorXd Tools::CalculateRMSE(const std::vector<VectorXd> &estimations, const std::vector<VectorXd> &ground_truth) {
	/**
	TODO:
	 * Calculate the RMSE here.
	 */
	VectorXd rmse(4);
	rmse.fill(0);

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 1; i < estimations.size(); ++i) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

double Tools::CalculeNIS(const VectorXd &Dz, const MatrixXd &Si) {
	//hold nis error
	double epislon = 0;
	//compute error
	epislon = Dz.transpose() * Si * Dz;
	//return error
	return epislon;
}

void Tools::PrintNIS(MeasurementPackage::SensorType sensor, std::vector<double>& values) {
	//common string values
	const std::string xTitle = "time";
	const std::string yTitle = "NIS values";

	std::string title;
	std::string fileName;

	if (sensor == MeasurementPackage::RADAR) {
		//print RADAR nis values
		title = "RADAR NIS Values";
		fileName = "nis_radar.png";
	} else {
		title = "LIDAR NIS Values";
		fileName = "nis_lidar.png";
	}
	Tools::PrintGraph(values, fileName, 5.991, title, xTitle, yTitle);
}

void Tools::PrintGraph(const std::vector<double> &data, const std::string fileName, const double reference,
		const std::string title, const std::string xTitle, const std::string yTitle) {

	//open a pipe to gnuplot
	static FILE *gnuplotPipe = NULL;
	if (gnuplotPipe == NULL) {
		gnuplotPipe = popen("gnuplot -persist", "w");
	}

	//print graph
	if (gnuplotPipe) {
		//prepare graph
		fprintf(gnuplotPipe, "reset\n"); //gnuplot commands    
		fprintf(gnuplotPipe, "set term png #output terminal and file\n");
		fprintf(gnuplotPipe, "set output '%s'\n", fileName.c_str());
		//fprintf(gnuplotPipe, "set xrange [min:max]\n");
		//fprintf(gnuplotPipe, "set yrange [0:]\n");
		fprintf(gnuplotPipe, "set style fill solid 0.5\n");
		fprintf(gnuplotPipe, "set xlabel '%s'\n", xTitle.c_str());
		fprintf(gnuplotPipe, "set ylabel '%s'\n", yTitle.c_str());
		fprintf(gnuplotPipe, "set title '%s'\n", title.c_str());
		//fprintf(gnuplotPipe, "plot(x, sin(x))\n");
		fprintf(gnuplotPipe, "plot '-' using 2:1 title 'NIS value' with linespoint, %f  title '95% Reference'\n", reference);
		unsigned int count = 0;
		for (double nis : data) {
			//ignore first point 
			if (count < 1) {
				count++;
				continue;
			}
			fprintf(gnuplotPipe, "%f ", nis);
			fprintf(gnuplotPipe, "%d ", count);
			fprintf(gnuplotPipe, "\n");
			count++;
		}

		//flush pipe
		fflush(gnuplotPipe);

		//close pipe
		fprintf(gnuplotPipe, "exit \n");
		pclose(gnuplotPipe);
		gnuplotPipe = NULL;
		std::cout << "Grap Plot finished\n";
	}
}