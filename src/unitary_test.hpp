#ifndef SRC_UNITARY_TEST_HPP_
#define SRC_UNITARY_TEST_HPP_

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
#include "Cartesian2Spherical.hpp"
#include "Estimation.hpp"

void triplet_cartesian (const std::vector<std::vector<double>> &t_match,
						const double &t_width,
						const double &t_height,
						Vec_Points<double> &sph) {


	for (auto &x:t_match) {
		Points<double> p = convertCartesian2Spherical (x[0], x[1], t_width, t_height);
		sph.push_back(p);
	}

}

void unitary_check(std::string directory_name, double err_tol) {

	std::string file_name = directory_name + "/data/20181218-160912-843294_20181218-160924-593294_20181218-160929-593290";

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
	triplet_cartesian(cart1, width, height, sph1);
	triplet_cartesian(cart2, width, height, sph2);
	triplet_cartesian(cart3, width, height, sph3);

	// construct the vector of sphere coordinates
	spheres.push_back(sph1);
	spheres.push_back(sph2);
	spheres.push_back(sph3);

	Vec_Points<double> sv_scene { };
	std::vector<Points<double>> positions { };

	std::vector<Mat_33<double>> sv_r_liste {};
	std::vector<Points<double>> sv_t_liste {};

	pose_estimation(spheres, err_tol, sv_scene, positions, sv_r_liste, sv_t_liste);

	// setup path to save the output result
	std::string output_file = directory_name + "/data/sv_scene.m";

	if (sv_scene.save_vecpoints(output_file)) {
		// Error opening the file
		throw std::runtime_error("Error writing the point cloud.");
	}

}

#endif /* SRC_UNITARY_TEST_HPP_ */
