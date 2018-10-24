#ifndef SRC_ESTIMATION_HPP_
#define SRC_ESTIMATION_HPP_

#include <stdio.h>
#include <iostream>
#include <vector>

// Include files to use OpenCV API
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <experimental/filesystem>
#include <limits>
#include <unistd.h>
#include "Points.hpp"
#include "Mat_33.hpp"
#include "Vec_Points.hpp"

template <typename T>
void estimation_rot_trans (const std::vector<Vec_Points<T>> &p3d_liste, const std::vector<std::vector<T>> sv_u_liste,
						   std::vector<Mat_33<T>> &sv_r_liste, std::vector<Points<T>> &sv_t_liste)
// Takes as inputs p3d_liste and sv_u_liste,
// and generates outputs sv_r_liste, sv_t_liste
{

	int nb_sph { p3d_liste.size() };
	int nb_pts { p3d_liste[0].size() };

	std::vector<Vec_Points<T>> p3d_liste_exp { };

	for (int i {0}; i < nb_sph; ++i) {
		p3d_liste_exp.push_back (p3d_liste[i] * sv_u_liste[i]);
	}

	std::vector<Points<T>> sv_cent_liste { };
	for (int i {0}; i < nb_sph; ++i) {
		sv_cent_liste.push_back (p3d_liste_exp[i].mean());
	}

	// Calculates the distances to the centers of the vector of points
	std::vector<Vec_Points<T>> sv_diff_liste { };
	for (int i {0}; i < nb_sph; ++i) {
		sv_diff_liste.push_back (p3d_liste_exp[i] - sv_cent_liste[i]);
	}

	// Multiply vector of points transposed by vector of points
	std::vector<Mat_33<T>> sv_corr_liste {};
	for (int i {0}; i < (nb_sph - 1); ++i) {
		sv_corr_liste.push_back = sv_diff_liste[i] * sv_diff_liste[i+1];
	}

	// check size of sv_r_liste. If it is empty fill with zero Mat_33
	if (sv_r_liste.size() < (nb_sph - 1)) {
		Mat_33<T> zm {};
		for (int i{sv_r_liste.size()}; i < (nb_sph - 1); ++i) {
			sv_r_liste.push_back (zm);
		}
	}

	// check size of sv_t_liste. If it is empty fill with zero Points
	if (sv_t_liste.size() < (nb_sph -1)) {
		Points<T> zp {};
		for (int i{sv_t_liste.size()}; i < (nb_sph - 1); ++i) {
			sv_t_liste.push_back (zp);
		}
	}

	for (int i{0}; i < (nb_sph -1); ++i) {
		// Matrices for SVD computation
		Mat_33<T> svd_Ut{};
		Mat_33<T> svd_V{};
		Mat_33<T> sv_r{};

		sv_corr_liste[i].svd(svd_Ut, svd_V);
		sv_r_liste[i].svd_rotation(svd_V, svd_Ut);
		sv_t_liste[i] = sv_cent_liste[i+1] - (sv_r_liste[i] * sv_cent_liste[i]);

	}

}

template <typename T>
inline Points<T> intersection (const Mat_33<T> &c, const Mat_33<T> &azim){
// Takes as input matrices c and azim and returns a point

	// For accumulating the results
	std::vector<Mat_33<T>> v{};
	std::vector<Points<T>> vp{};

	for (int i{0}; i<3; ++i) {
		// takes row i from azim
		T a0 = azim[i][0];
		T a1 = azim[i][1];
		T a2 = azim[i][2];

		// computes identity(3) - azim[i]' * azim[i]
		Mat_33<T> v1{1 - a0*a0 ,    -a0*a1 ,    -a0*a2,
						-a1*a0 , 1 - a1*a1 ,    -a1*a2,
						-a2*a0 ,    -a2*a1 , 1 - a2*a2};

		v.push_back(v1);

		// takes the row of c
		Points<T> row {c[i][0], c[i][1], c[i][2]};

		Points<T> vp1 {v1 * row};

		vp.push_back(vp1);
	}

	Mat_33<T> sum_v = v[0] + v[1] + v[2];
	Points<T> sum_vp = vp[0] + vp[1] + vp[2];

	// Computes the inverse of the matrix
	Mat_33<T> sum_v_inv {sum_v.inv()};

	Points<T> inter {sum_v_inv * sum_vp};

	return inter;
}

template <typename T>
void centers_determination (const std::vector<Mat_33<T>> &sv_r_liste, const std::vector<Points<T>> &sv_t_liste,
							std::vector<Points<T>> &center_liste)
// Takes sv_r_liste and sv_t_liste as inputs
// Generates center_liste as output
{
	int nb_sph { sv_r_liste.size() + 1 };

	// if center_liste is empty
	if (center_liste.size() < nb_sph) {
		Points<T> zp{};
		for (int i{center_liste.size()}; i < nb_sph; ++i) {
			center_liste.push_back (zp);
		}
	}

	for (int i{0}; i < (nb_sph -1); ++i) {
		Points<T> center = {};
		for (int j{0}; j < (i-1); ++j) {
			int k { i - 1 - j };
			center = (sv_r_liste[k].transpose()) * (center - sv_t_liste[k]);
		}
		center_liste[i] = center;
	}

}

template <typename T>
std::vector<Points<T>> azim_determination (const std::vector<Points<T>> &azim_liste,
		                                   const std::vector<Mat_33<T>> &sv_r_liste,
										   const std::vector<Points<T>> &sv_t_liste)
{
	int nb_sph { azim_liste.size() };
	for (int i {0}; i < nb_sph; ++i) {
		for (int j{0}; j < i; ++j) {
			int k { i - 1 - j };
			azim_liste[i] = sv_r_liste[k].transpose() * azim_liste[i];
		}
	}
}

template<typename T>
void estimation_rayons(const std::vector<Vec_Points<T>> &p3d_liste, std::vector<std::vector<T>> sv_u_liste, const std::vector<Mat_33<T>> &sv_r_liste,
		const std::vector<Points<T>> &sv_t_liste, std::vector<T> &sv_e_liste) {
// Takes as input p3d_liste, sv_u_liste, sv_r_liste and sv_t_liste
// and generates as output sv_u_liste and sv_e_liste

	int nb_sph { p3d_liste.size() };
	int nb_pts { p3d_liste[0].size() };

	std::vector<Points<T>> center_liste { };
	centers_determination(sv_r_liste, sv_t_liste, center_liste);

	if (sv_e_liste.size() < (nb_sph - 1)) {
		for (int i { sv_e_liste.size() }; i < (nb_sph - 1); ++i) {
			sv_e_liste.push_back(0);
		}
	} else {
		for (int i { 0 }; i < (nb_sph - 1); ++i) {
			sv_e_liste[i] = 0;
		}
	}

	for (int j { 0 }; j < nb_pts; ++j) {
		std::vector<Points<T>> azim_liste { };
		for (int i { 0 }; i < nb_sph; ++i) {
			azim_liste.push_back(p3d_liste[i][j]);
		}
		azim_liste = azim_determination(azim_liste, sv_r_liste, sv_t_liste);
		std::vector<T> rayons { };
		try {
			rayons = intersection(center_liste, azim_liste);
			for (int k { 0 }; k < nb_sph; ++k) {
				sv_u_liste[k][j] = rayons[k];
			}
		} catch (...) {
			for (int k { 0 }; k < nb_sph; ++k) {
				rayons[k] = sv_u_liste[k][j];
			}
		}
		std::vector<Points<T>> inter_liste { };

		for (int i { 0 }; i < nb_sph; ++i) {
			inter_liste[i] = center_liste[i] + azim_liste[i] * rayons[i];
		}

		for (int i { 0 }; i < (nb_sph - 1); ++i) {
			Points<T> p { inter_liste[i] - inter_liste[i + 1] };
			T n { p.norm() };
			sv_e_liste[i] = sv_e_liste[i] > n ? sv_e_liste[i] : n;
		}

	}
}

template <typename T>
void pose_scene (const Vec_Points<T> &p3d_1, const Vec_Points<T> &p3d_2, const Vec_Points<T> &p3d_3,
				 const Mat_33<T> &sv_r_12, const Mat_33<T> &sv_r_23, const Mat_33<T> &sv_r_31,
				 const Points<T> &sv_t_12, const Points<T> &sv_t_23, const Points<T> &sv_t_31,
				 Vec_Points<T> &sv_scene) {

	size_t longueur {p3d_1.size()};

	Points<T> c1 {0, 0, 0};
	Points<T> c2 {sv_t_12};	// c2 = c1 + sv_t_12
	Points<T> c3 {sv_t_12 + sv_r_12 * sv_t_23};	// c3 = c2 + sv_r_12 * sv_t_23

	Vec_Points<T> azim1m {p3d_1};
	Vec_Points<T> azim2m {(p3d_2 * sv_r_23) * sv_r_31};
	Vec_Points<T> azim3m {p3d_3 * sv_r_31};

	for (size_t i{0}; i < longueur; ++i) {

		Points<T> azim1 {azim1m[i]};
		Points<T> azim2 {azim2m[i]};
		Points<T> azim3 {azim3m[i]};

		Mat_33<T> c {c1, c2, c3};
		Mat_33<T> azim { azim1, azim2, azim3 };

		Points<T> inter{};
		inter = intersection (c, azim);

		sv_scene[i] = inter;
	}

}

template <typename T>
void pose_estimation (const std::vector<Vec_Points<double>> &p3d_liste, const double error_max,
					  Vec_Points<T> &sv_scene,
					  Points<double> &positions) {

// modified for n-tuple

	int nb_sph = p3d_liste.size();
	int nb_pts = p3d_liste[0].size();

	std::vector<std::vector<T>> sv_u_liste {};
	std::vector<Points<T>> sv_t_liste {};
	std::vector<Mat_33<T>> sv_r_liste {};
	std::vector<T> sv_e_liste {};

	std::vector<T> ones(nb_pts,1);
	for (int i=0; i < nb_sph; ++i) {
		sv_u_liste.push_back(ones);
	}

	T sv_e_old { 0 };
	T sv_e_norm { 1 };
	int count { 0 };

	while (abs(sv_e_norm - sv_e_old) > error_max) {

		estimation_rot_trans (p3d_liste, sv_u_liste,
							  sv_r_liste, sv_t_liste);

		estimation_rayons (p3d_liste, sv_u_liste, sv_r_liste, sv_t_liste,
						   sv_u_liste, sv_e_liste);

		++count;
		T sv_t_norm { 0 };
		for (int i {0}; i < sv_t_liste.size(); ++i) {
			sv_t_norm += (sv_t_liste[i].norm());
		}

		T max { 0 };
		for (const auto x: sv_e_liste) {
			max = x > max ? x : max;
		}
		sv_e_norm = sv_e_liste.size() * max / sv_t_norm;
	}

	pose_scene (p3d_liste, sv_u_liste, sv_r_liste, sv_t_liste,
			    sv_scene, positions);

}

#endif /* SRC_ESTIMATION_HPP_ */
