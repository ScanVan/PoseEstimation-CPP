#ifndef SRC_TEST_TRIPLET_HPP_
#define SRC_TEST_TRIPLET_HPP_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <exception>

#include "Vec_Points.hpp"
#include <opencv4/opencv2/core/types.hpp>
#include <opencv4/opencv2/core.hpp>
#include <opencv4/opencv2/features2d.hpp>
#include <opencv4/opencv2/imgcodecs.hpp>
#include <opencv4/opencv2/opencv.hpp>
#include <math.h>

template <typename T>
void print1(const std::vector<T> &v) {
	for (auto const &x: v ) {
		std::cout << std::setprecision(15) <<  x << " ";
	}
	std::cout << std::endl;
}

void triplet_cartesian (const std::vector<std::vector<double>> &t_match,
						const double &t_width,
						const double &t_height,
						Vec_Points<double> &sph) {


	for (auto &x:t_match) {
		// coordinate re-normalization
		double tm1 = (x[0] / t_width) * 2 * M_PI;
		double tm2 = (0.5 - (x[1] / t_height)) * M_PI;

		// coordinate conversion
		double p1 { cos (tm2) * cos (tm1) };
		double p2 { cos (tm2) * sin (tm1) };
		double p3 { sin(tm2) };

		Points<double> p {p1, p2, p3};
		sph.push_back(p);
	}

}


void first_check(std::string directory_name, double err_tol) {

	std::string file_name = directory_name + "/data/triplet_matches";

	std::ifstream file(file_name);

	// check if file is open
	if (!file.good()) {
		throw (std::runtime_error ("Input file could not be opened!"));
	}

	// create the spheres
	std::vector<Vec_Points<double>> spheres { };
	Vec_Points<double> sph1 { };
	Vec_Points<double> sph2 { };
	Vec_Points<double> sph3 { };

	// vector of cartesian coordinates
	std::vector<std::vector<double>> cart1{};
	std::vector<std::vector<double>> cart2{};
	std::vector<std::vector<double>> cart3{};

	// read from file the cartesian coordinates
	while (file.good()) {
		double x1 { }, y1 { }, x2 { }, y2 { }, x3 { }, y3 { };
		file >> std::scientific >> std::setprecision(15) >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
		if (file.good()) {
			std::vector<double> p1 { x1, y1 };
			cart1.push_back(p1);
			std::vector<double> p2 { x2, y2 };
			cart2.push_back(p2);
			std::vector<double> p3 { x3, y3 };
			cart3.push_back(p3);
		}
	}

	double width = 6016;
	double height = 3008;

	// convert to unitary sphere coordinates
	triplet_cartesian (cart1, width, height, sph1);
	triplet_cartesian (cart2, width, height, sph2);
	triplet_cartesian (cart3, width, height, sph3);

	// construct the vector of sphere coordinates
	spheres.push_back(sph1);
	spheres.push_back(sph2);
	spheres.push_back(sph3);

	Vec_Points<double> sv_scene { };
	std::vector<Points<double>> positions { };

	pose_estimation (spheres, err_tol, sv_scene, positions);

	// setup path to save the output result
	std::string output_file = directory_name +  "/data/sv_scene.m";

	if (sv_scene.save_vecpoints(output_file)) {
		// Error opening the file
		throw std::runtime_error("Error writing the point cloud.");
	}

}



#endif /* SRC_TEST_TRIPLET_HPP_ */
