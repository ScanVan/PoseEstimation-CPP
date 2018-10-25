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
#include "Print_Data.hpp"
#include <cmath>

template <typename T>
void estimation_rot_trans (const std::vector<Vec_Points<T>> &p3d_liste, const std::vector<std::vector<T>> sv_u_liste,
						   std::vector<Mat_33<T>> &sv_r_liste, std::vector<Points<T>> &sv_t_liste)
// Takes as inputs p3d_liste and sv_u_liste,
// and generates outputs sv_r_liste, sv_t_liste
{

	size_t nb_sph { p3d_liste.size() };
	//size_t nb_pts { p3d_liste[0].size() };

	std::vector<Vec_Points<T>> p3d_liste_exp { };

	for (size_t i { 0 }; i < nb_sph; ++i) {
		p3d_liste_exp.push_back(p3d_liste[i] * sv_u_liste[i]);
	}

	std::vector<Points<T>> sv_cent_liste { };
	for (size_t i { 0 }; i < nb_sph; ++i) {
		sv_cent_liste.push_back(p3d_liste_exp[i].mean());
	}

	// Calculates the distances to the centers of the vector of points
	std::vector<Vec_Points<T>> sv_diff_liste { };
	for (size_t i { 0 }; i < nb_sph; ++i) {
		sv_diff_liste.push_back(p3d_liste_exp[i] - sv_cent_liste[i]);
	}

	// Multiply vector of points transposed by vector of points
	std::vector<Mat_33<T>> sv_corr_liste { };
	for (size_t i { 0 }; i < (nb_sph - 1); ++i) {
		sv_corr_liste.push_back(sv_diff_liste[i] * sv_diff_liste[i + 1]);
	}

	// check size of sv_r_liste. If it is empty fill with zero Mat_33
	if (sv_r_liste.size() < (nb_sph - 1)) {
		Mat_33<T> zm { };
		for (size_t i { sv_r_liste.size() }; i < (nb_sph - 1); ++i) {
			sv_r_liste.push_back(zm);
		}
	}

	// check size of sv_t_liste. If it is empty fill with zero Points
	if (sv_t_liste.size() < (nb_sph - 1)) {
		Points<T> zp { };
		for (size_t i { sv_t_liste.size() }; i < (nb_sph - 1); ++i) {
			sv_t_liste.push_back(zp);
		}
	}

	for (size_t i { 0 }; i < (nb_sph - 1); ++i) {
		// Matrices for SVD computation
		Mat_33<T> svd_Ut { };
		Mat_33<T> svd_V { };
		Mat_33<T> sv_r { };

		sv_corr_liste[i].svd(svd_Ut, svd_V);
		sv_r_liste[i].svd_rotation(svd_V, svd_Ut);
		sv_t_liste[i] = sv_cent_liste[i + 1] - (sv_r_liste[i] * sv_cent_liste[i]);

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
std::vector<T> intersection_bis (const std::vector<Points<T>> &liste_p, const std::vector<Points<T>> &liste_azim ){

	size_t nb_pts { liste_p.size() };
	Mat_33<T> sum_v { };
	Points<T> sum_vp { };

	for (size_t i { 0 }; i < nb_pts; ++i) {
		T v11 { 1.0 - liste_azim[i][0] * liste_azim[i][0] };
		T v22 { 1.0 - liste_azim[i][1] * liste_azim[i][1] };
		T v33 { 1.0 - liste_azim[i][2] * liste_azim[i][2] };
		T v12 { -liste_azim[i][0] * liste_azim[i][1] };
		T v13 { -liste_azim[i][0] * liste_azim[i][2] };
		T v23 { -liste_azim[i][1] * liste_azim[i][2] };
		sum_v[0][0] = sum_v[0][0] + v11;
		sum_v[0][1] += v12;
		sum_v[0][2] += v13;
		sum_v[1][0] += v12;
		sum_v[1][1] += v22;
		sum_v[1][2] += v23;
		sum_v[2][0] += v13;
		sum_v[2][1] += v23;
		sum_v[2][2] += v33;
		T p1 { liste_p[i][0] };
		T p2 { liste_p[i][1] };
		T p3 { liste_p[i][2] };
		sum_vp[0] += p1 * v11 + p2 * v12 + p3 * v13;
		sum_vp[1] += p1 * v12 + p2 * v22 + p3 * v23;
		sum_vp[2] += p1 * v13 + p2 * v23 + p3 * v33;
	}

	Points<T> inter { sum_v.inv() * sum_vp };

	std::vector<T> rayons { };
	for (size_t i { 0 }; i < nb_pts; ++i) {
		rayons.push_back(0);
	}

	for (size_t i { 0 }; i < nb_pts; ++i) {
		Points<T> centre { liste_p[i] };
		Points<T> azim { liste_azim[i] };
		Points<T> inter_proj { azim * ((inter - centre) * azim) / (azim * azim) };
		T direction { inter_proj * azim };
		if (direction < 0) {
			rayons[i] = -inter_proj.norm();
		} else {
			rayons[i] = inter_proj.norm();
		}
	}
	return rayons;
}


template <typename T>
void centers_determination (const std::vector<Mat_33<T>> &sv_r_liste, const std::vector<Points<T>> &sv_t_liste,
							std::vector<Points<T>> &center_liste)
// Takes sv_r_liste and sv_t_liste as inputs
// Generates center_liste as output
{
	size_t nb_sph { sv_r_liste.size() + 1 };

	// if center_liste is empty
	if (center_liste.size() < nb_sph) {
		Points<T> zp { };
		for (size_t i { center_liste.size() }; i < nb_sph; ++i) {
			center_liste.push_back(zp);
		}
	}

	for (size_t i { 0 }; i < nb_sph; ++i) {
		Points<T> center = { };
		for (int j { 0 }; j < static_cast<int>(i); ++j) {
			size_t k { i - 1 - j };
			center = (sv_r_liste[k].transpose()) * (center - sv_t_liste[k]);

		}
		center_liste[i] = center;
	}
}

template<typename T>
std::vector<Points<T>> azim_determination(std::vector<Points<T>> &azim_liste,
										  const std::vector<Mat_33<T>> &sv_r_liste,
										  const std::vector<Points<T>> &sv_t_liste) {

	size_t nb_sph { azim_liste.size() };
	for (size_t i { 0 }; i < nb_sph; ++i) {
		for (size_t j { 0 }; j < i; ++j) {
			size_t k { i - 1 - j };
			azim_liste[i] = sv_r_liste[k].transpose() * azim_liste[i];
		}
	}

	return azim_liste;

}


template<typename T>
void estimation_rayons(const std::vector<Vec_Points<T>> &p3d_liste, std::vector<std::vector<T>> &sv_u_liste, const std::vector<Mat_33<T>> &sv_r_liste,
		const std::vector<Points<T>> &sv_t_liste, std::vector<T> &sv_e_liste) {
// Takes as input p3d_liste, sv_u_liste, sv_r_liste and sv_t_liste
// and generates as output sv_u_liste and sv_e_liste

	size_t nb_sph { p3d_liste.size() };
	size_t nb_pts { p3d_liste[0].size() };

	std::vector<Points<T>> center_liste { };

	// Initialize center_liste
	if (center_liste.size() < nb_sph) {
		Points<T> p {};
		for (size_t i{ center_liste.size() }; i < nb_sph; ++i ) {
			center_liste.push_back(p);
		}
	}

	centers_determination(sv_r_liste, sv_t_liste, center_liste);

	// Initialize sv_e_liste
	if (sv_e_liste.size() < (nb_sph - 1)) {
		for (size_t i { sv_e_liste.size() }; i < (nb_sph - 1); ++i) {
			sv_e_liste.push_back(0);
		}
	} else {
		for (size_t i { 0 }; i < (nb_sph - 1); ++i) {
			sv_e_liste[i] = 0;
		}
	}

	for (size_t j { 0 }; j < nb_pts; ++j) {

		std::vector<Points<T>> azim_liste { };
		for (size_t i { 0 }; i < nb_sph; ++i) {
			azim_liste.push_back(p3d_liste[i][j]);
		}

		azim_determination(azim_liste, sv_r_liste, sv_t_liste);

		std::vector<T> rayons { };
		try {
			rayons = intersection_bis(center_liste, azim_liste);

			for (size_t k { 0 }; k < nb_sph; ++k) {
				sv_u_liste[k][j] = rayons[k];
			}
		} catch (...) {
			for (size_t k { 0 }; k < nb_sph; ++k) {
				rayons[k] = sv_u_liste[k][j];
			}
		}

		std::vector<Points<T>> inter_liste { };

		// Initialize inter_liste
		for (size_t i { 0 }; i < nb_sph; ++i) {
			Points<T> p {};
			inter_liste.push_back(p);
		}

		for (size_t i { 0 }; i < nb_sph; ++i) {
			inter_liste[i] = center_liste[i] + azim_liste[i] * rayons[i];
		}

		for (size_t i { 0 }; i < (nb_sph - 1); ++i) {
			Points<T> p { inter_liste[i] - inter_liste[i + 1] };
			T n { p.norm() };
			sv_e_liste[i] = sv_e_liste[i] > n ? sv_e_liste[i] : n;
		}

	}
}

template <typename T>
void pose_scene (const std::vector<Vec_Points<T>> &p3d_liste,
			     const std::vector<std::vector<T>> &sv_u_liste,
			     const std::vector<Mat_33<T>> &sv_r_liste,
				 const std::vector<Points<T>> &sv_t_liste,
				 Vec_Points<T> &sv_scene,
				 std::vector<Points<T>> &center_liste) {

	size_t nb_sph { p3d_liste.size() };
	size_t nb_pts { p3d_liste[0].size() };

	centers_determination(sv_r_liste, sv_t_liste, center_liste);

	// If sv_scene empty, fill with zeros
	if (sv_scene.size() < nb_pts) {
		for (size_t i{0}; i < nb_pts; ++i) {
			sv_scene.push_back(0, 0, 0);
		}
	}

	for (size_t j { 0 }; j < nb_pts; ++j) {
		std::vector<Points<T>> azim_liste { };
		for (size_t i { 0 }; i < nb_sph; ++i) {
			azim_liste.push_back(p3d_liste[i][j]);
		}
		azim_liste = azim_determination(azim_liste, sv_r_liste, sv_t_liste);
		std::vector<T> rayons { };
		try {
			rayons = intersection_bis(center_liste, azim_liste);
		} catch (...) {
			for (size_t k { 0 }; k < nb_sph; ++k) {
				rayons[k] = sv_u_liste[k][j];
			}
		}
		std::vector<Points<T>> inter_liste { };

		// Initialize inter_liste
		if (inter_liste.size() < nb_sph) {
			for (size_t i { inter_liste.size() }; i < nb_sph; ++i) {
				Points<T> p { };
				inter_liste.push_back(p);
			}
		}

		for (size_t i { 0 }; i < nb_sph; ++i) {
			inter_liste[i] = center_liste[i] + azim_liste[i] * rayons[i];
		}

		Points<T> inter { };
		for (size_t i {0}; i < inter_liste.size(); ++i) {
			inter = inter + inter_liste[i];
		}
		sv_scene[j] = inter / inter_liste.size();

	}






	/*size_t longueur {p3d_1.size()};

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
	}*/

}


template <typename T>
void pose_estimation (const std::vector<Vec_Points<T>> &p3d_liste, const T error_max,
					  Vec_Points<T> &sv_scene,
					  std::vector<Points<T>> &positions) {

// modified for n-tuple

	size_t nb_sph = p3d_liste.size();
	size_t nb_pts = p3d_liste[0].size();

	std::vector<std::vector<T>> sv_u_liste {};
	std::vector<Points<T>> sv_t_liste {};
	std::vector<Mat_33<T>> sv_r_liste {};
	std::vector<T> sv_e_liste {};

	// Initialize sv_u_liste
	std::vector<T> ones(nb_pts, 1);
	for (size_t i { 0 }; i < nb_sph; ++i) {
		sv_u_liste.push_back(ones);
	}

	// Initialize sv_r_liste
	for (size_t i { 0 }; i < (nb_sph - 1); ++i) {
		Mat_33<T> m { };
		sv_r_liste.push_back(m);
	}

	// Initialize sv_t_liste;
	for (size_t i { 0 }; i < (nb_sph - 1); ++i) {
		Points<T> p { };
		sv_t_liste.push_back(p);
	}

	T sv_e_old { 0 };
	T sv_e_norm { 1 };
	int count { 0 };

	T diff_error { sv_e_norm - sv_e_old };

	while (diff_error > error_max) {

		sv_e_old = sv_e_norm;
		count ++;

		estimation_rot_trans (p3d_liste, sv_u_liste,
							  sv_r_liste, sv_t_liste);

		estimation_rayons (p3d_liste, sv_u_liste, sv_r_liste, sv_t_liste, sv_e_liste);

		T sv_t_norm { 0 };
		for (size_t i {0}; i < sv_t_liste.size(); ++i) {
			sv_t_norm += (sv_t_liste[i].norm());
		}

		T max_num { sv_e_liste[0] };
		for (size_t i { 1 }; i < sv_e_liste.size(); ++i) {
			if (sv_e_liste[i] > max_num)
				max_num = sv_e_liste[i];
		}

		sv_e_norm = sv_e_liste.size() * max_num / sv_t_norm;
		diff_error = (sv_e_norm - sv_e_old) > 0 ? (sv_e_norm - sv_e_old):(sv_e_old - sv_e_norm);

	}

	// Initialize positions
	if (positions.size() < nb_sph) {
		for (size_t i { 0 }; i < nb_sph; ++i) {
			Points<T> p {};
			positions.push_back(p);
		}
	}

	pose_scene (p3d_liste, sv_u_liste, sv_r_liste, sv_t_liste,
			    sv_scene, positions);

}

#endif /* SRC_ESTIMATION_HPP_ */
