//============================================================================
// Name        : PoseEstimation.cpp
// Author      : Marcelo Kaihara
// Version     : 2.0
// Copyright   : 
// Description : Pose Estimation algorithm implemented in C++
//============================================================================

#include <iostream>

#include "Mat_33.hpp"
#include "Points.hpp"
#include "Vec_Points.hpp"
#include "Estimation.hpp"
#include <chrono>
#include "Print_Data.hpp"
//#include "test_prior_model_nuple.hpp"
#include "test_triplet.hpp"

std::string GetCurrentWorkingDir( void ) {
// gets the current working directory
  char buff[FILENAME_MAX]{};
  if (getcwd( buff, FILENAME_MAX )==nullptr) {
	  throw (std::runtime_error("The directory could not be determined."));
  }
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main() {
	std::string file_name {GetCurrentWorkingDir() +"/data/triplet_matches"};
	first_check(file_name, 1e-8);

}

/*
int main() {

	std::cout << "simulation(42,3,1e-8);" << std::endl;
	simulation(42,3,1e-8);
	std::cout << "simulation(42,6,1e-8);" << std::endl;
	simulation(42,6,1e-8);
	std::cout << "simulation(42,9,1e-8);" << std::endl;
	simulation(42,9,1e-8);


	std::cout << "simulation(16,3,1e-8);" << std::endl;
	simulation(16,3,1e-8);
	std::cout << "simulation(16,6,1e-8);" << std::endl;
	simulation(16,6,1e-8);
	std::cout << "simulation(16,9,1e-8);" << std::endl;
	simulation(16,9,1e-8);

	std::cout << "simulation(19,3,1e-8);" << std::endl;
	simulation(10,3,1e-8);
	std::cout << "simulation(10,6,1e-8);" << std::endl;
	simulation(10,6,1e-8);
	std::cout << "simulation(10,9,1e-8);" << std::endl;
	simulation(10,9,1e-8);

	std::cout << "end" << std::endl;

	return 0;

}
*/

/*
int main2() {

	// Set path to input data
	std::string path2data = GetCurrentWorkingDir() + "/data/";
	std::string path_data1 = path2data + "p3d_1.txt";
	std::string path_data2 = path2data + "p3d_2.txt";
	std::string path_data3 = path2data + "p3d_3.txt";

	// Input vector of points
	Vec_Points<double> p3d_1 { };
	Vec_Points<double> p3d_2 { };
	Vec_Points<double> p3d_3 { };

	// Load data from file
	if (p3d_1.load_vecpoints(path_data1)) {
		// Error opening the file
		return 1;
	}
	if (p3d_2.load_vecpoints(path_data2)) {
		// Error opening the file
		return 1;
	}
	if (p3d_3.load_vecpoints(path_data3)) {
		// Error opening the file
		return 1;
	}

	// Output result, vector of points
	Vec_Points<double> sv_scene{p3d_1.size()};

	// Output result, position of the center
	std::vector<Points<double>> positions {};

	// Input. List of vec points
	std::vector<Vec_Points<double>> p3d_liste { p3d_1, p3d_2, p3d_3 };

	// Input. error
	double error_max { 1e-8 };

	// Counts the number of iterations
	size_t num_iter { 0 };

	// For timing measurements
	std::chrono::high_resolution_clock::time_point t1{};
	std::chrono::high_resolution_clock::time_point t2{};

	// start measuring time
	t1 = std::chrono::high_resolution_clock::now();

	// main algorithm
	pose_estimation (p3d_liste, error_max,
					 sv_scene, positions, num_iter);

	// stop measuring time
	t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	std::cout << "Number of points: " << p3d_1.size() << std::endl;
	std::cout << "Number of iterations: " << num_iter << std::endl;
	std::cout << "Execution time: " << duration << " microseconds" << std::endl;

	// setup path to save the output result
	std::string path_out_data1 = path2data + "sv_scene_cpp.m";

	if (sv_scene.save_vecpoints(path_out_data1)) {
		// Error opening the file
		return 1;
	}



	return 0;
}
*/
