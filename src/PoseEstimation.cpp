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
#include "unitary_test.hpp"

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
	std::string directory_name {GetCurrentWorkingDir()};
	unitary_check (directory_name, 1e-8);
}
