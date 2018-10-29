#ifndef SRC_TEST_PRIOR_MODEL_NUPLE_HPP_
#define SRC_TEST_PRIOR_MODEL_NUPLE_HPP_

#include "Vec_Points.hpp"
#include <opencv4/opencv2/core/types.hpp>
#include <opencv4/opencv2/core.hpp>
#include <opencv4/opencv2/features2d.hpp>
#include <opencv4/opencv2/imgcodecs.hpp>
#include <opencv4/opencv2/opencv.hpp>


void importation_data(std::string path, Vec_Points<double> &data) {

	std::ifstream ifs (path, std::ifstream::in);

	if (!ifs.is_open())
		throw std::runtime_error ("Data file cannot be opened.");

	std::string line { };
	while (getline(ifs, line)) {
		double d1 { }, d2 { }, d3 { };
		std::stringstream ss { };
		ss << line;
		ss >> d1 >> d2 >> d3;
		data.push_back(d1, d2, d3);
	}
}

void importation_centers (std::string path, int num, Vec_Points<double> &data) {

	std::ifstream ifs (path, std::ifstream::in);

	if (!ifs.is_open())
		throw std::runtime_error ("Data file with centers cannot be opened.");

	std::string line { };
	int counter { 0 };
	while (getline(ifs, line) && (counter < num)) {
		double d1 { }, d2 { }, d3 { };
		std::stringstream ss { };
		ss << line;
		ss >> d1 >> d2 >> d3;
		data.push_back(d1, d2, d3);
		counter++;
	}

}

void projection (Vec_Points<double> &p3d, Points<double> &center, Vec_Points<double> &res) {

	for (size_t i { 0 }; i < p3d.size(); ++i) {
		Points<double> vec { p3d[i] - center };
		vec = vec / (vec.norm());
		res.push_back(vec);
	}
}

void generate_unit_sphere (Vec_Points<double> &p3d, Vec_Points<double> &centers, std::vector<Vec_Points<double>> &spheres) {

	for (size_t i {0}; i < centers.size(); ++i) {
		Vec_Points<double> p3d_proj {};
		projection (p3d, centers[i], p3d_proj);
		spheres.push_back(p3d_proj);
	}
}


void writePly(std::string file, Vec_Points<double> features){
//	‘__gnu_cxx::__normal_iterator<ModelKeypoint*, std::vector<ModelKeypoint> >
	std::ofstream s {};
	s.open (file);

	s << "ply" << std::endl;
	s << "format ascii 1.0 " << std::endl;
	s << "element vertex " << features.size() << std::endl;
	s << "comment " << file << std::endl;
	s << "property float32 x " << std::endl;
	s << "property float32 y " << std::endl;
	s << "property float32 z " << std::endl;
	s << "end_header " << std::endl;

	for(size_t i {0}; i < features.size(); ++i){
		Points<double> f = features[i];
		s << f[0] << " " << f[1] << " " << f[2] << std::endl;
	}

	s.close();
}


void simulation (int num_model, int num_cam, double stop_crit) {
	std::string dirname {std::to_string(num_model)};
	std::string folder { "./prior_generated_test_models/" };
	Vec_Points<double> data {};
	importation_data (folder + dirname + "/model.dat", data);
	Vec_Points<double> centers_ori {};
	importation_centers (folder + dirname + "/centers_R1.txt", num_cam, centers_ori);
	std::vector<Vec_Points<double>> spheres {};
	generate_unit_sphere(data, centers_ori, spheres);
	Vec_Points<double> sv_scene {};
	std::vector<Points<double>> positions {};
	pose_estimation (spheres, stop_crit, sv_scene, positions);
	std::string plyFileName { "./models_est_" + dirname + "_" + std::to_string(num_cam) + ".ply"};
	writePly (plyFileName, sv_scene);
}



#endif /* SRC_TEST_PRIOR_MODEL_NUPLE_HPP_ */
