/*
 * util.hpp
 *
 *  Created on: 24/12/2013
 *      Author: polymorpher
 */

#pragma once
#include<vector>
#include<random>
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<cstring>
#include<ctime>
#include<unistd.h>
#include<cassert>
#include<google/dense_hash_map>
#include<memory>
#include<thread>
#include<sstream>
#define __DEBUG
#define CLOG(X) if(__DEBUG)printf(X)
typedef std::vector<std::vector<int> > Vector2DInt;
typedef std::vector<std::vector<double> > Vector2DDouble;
typedef std::vector<std::vector<long double> > Vector2DLongDouble;
typedef std::vector<int> VecInt;
typedef std::vector<double> VecDouble;
typedef std::vector<long double> VecLongDouble;
typedef google::dense_hash_map<int,int> MapIntInt;
typedef google::dense_hash_map<int,double> MapIntDouble;
typedef std::vector<MapIntInt> VecMapIntInt;
typedef std::vector<MapIntDouble> VecMapIntDouble;

inline void vecResize2D(Vector2DInt &a, int s0, int s1) {
	a.resize(s0);
	for (size_t i = 0; i < a.size(); i++) {
		a[i].resize(s1);
	}
}
inline void vecResize2D(Vector2DDouble &a, int s0, int s1) {
	a.resize(s0);
	for (size_t i = 0; i < a.size(); i++) {
		a[i].resize(s1);
	}
}

template<class T> T convertString(std::string s){
	std::stringstream ss;
	ss<<s;
	T t;
	ss>>t;
	return t;
}