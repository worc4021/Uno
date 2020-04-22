#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include "Logger.hpp"

std::vector<double> add_vectors(std::vector<double>& x, std::vector<double>& y, double scaling_factor);

double norm(std::vector<double>& x, double chosen_norm);

double norm_1(std::vector<double>& x);
double norm_1(std::map<int,double>& x);
double norm_1(std::vector<std::map<int, double> >& m);

double norm_2_squared(std::vector<double>& x);
double norm_2(std::vector<double>& x);

double norm_inf(std::vector<double>& x, unsigned int length = std::numeric_limits<unsigned int>::max());
double norm_inf(std::vector<std::map<int, double> >& m);

double dot(std::vector<double>& x, std::vector<double>& y);
double dot(std::vector<double>& x, std::map<int,double>& y);
double dot(std::map<int,double>& x, std::map<int,double>& y);

template <typename T>
void print_vector(std::ostream &stream, std::vector<T> x, unsigned int start = 0, unsigned int length = std::numeric_limits<unsigned int>::max()) {
    for (unsigned int i = start; i < std::min<unsigned int>(start + length, x.size()); i++) {
            stream << x[i] << " ";
    }
    stream << "\n";
    return;
}

template <typename T>
void print_vector(const Level& level, std::vector<T> x, unsigned int start = 0, unsigned int length = std::numeric_limits<unsigned int>::max()) {
    for (unsigned int i = 0; i < std::min<unsigned int>(start + length, x.size()); i++) {
            level << x[i] << " ";
    }
    level << "\n";
    return;
}

template <typename T, typename U>
void print_vector(std::ostream &stream, std::map<T, U> x) {
    for (std::pair<T, U> element: x) {
        T i = element.first;
        U xi = element.second;
        stream << "x[" << i << "] = " << xi << ", ";
    }
    stream << "\n";
    return;
}

template <typename T, typename U>
void print_vector(const Level& level, std::map<T, U> x) {
    for (std::pair<T, U> element: x) {
        T i = element.first;
        U xi = element.second;
        level << "x[" << i << "] = " << xi << ", ";
    }
    level << "\n";
    return;
}

#endif // UTILS_H
