#include <cmath>
#include "Utils.hpp"

std::vector<double> add_vectors(std::vector<double>& x, std::vector<double>& y, double scaling_factor) {
	if (x.size() != y.size()) {
		throw std::length_error("Utils.add_vectors: x and y have different sizes");
	}
	
	std::vector<double> z(x.size());
	for (unsigned int i = 0; i < x.size(); i++) {
		z[i] = x[i] + scaling_factor*y[i];
	}
	return z;
}

/* compute ||x||_1 */
double norm_1(std::vector<double>& x) {
	double norm = 0.;
	for (unsigned int i = 0; i < x.size(); i++) {
		norm += std::abs(x[i]);
	}
	return norm;
}

double norm_1(std::map<int,double>& x) {
	double norm = 0.;
	for (std::map<int,double>::iterator it = x.begin(); it != x.end(); it++) {
		double xi = it->second;
		norm += std::abs(xi);
	}
	return norm;
}

/* compute ||x||_2 */
double norm_2(std::vector<double>& x) {
	double norm = 0.;
	for (unsigned int i = 0; i < x.size(); i++) {
		norm += x[i]*x[i];
	}
	norm = std::sqrt(norm);
	return norm;
}

/* compute ||x||_infty */
double norm_inf(std::vector<double>& x) {
	double norm = 0.;
	for (unsigned int i = 0; i < x.size(); i++) {
		norm = std::max(norm, std::abs(x[i]));
	}
	return norm;
}

double dot(std::vector<double>& x, std::vector<double>& y) {
	if (x.size() != y.size()) {
		throw std::length_error("Utils.dot: x and y have different sizes");
	}
	
	double dot = 0.;
	for (unsigned int i = 0; i < x.size(); i++) {
		dot += x[i]*y[i];
	}
	return dot;
}