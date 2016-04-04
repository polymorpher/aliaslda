#ifndef GeneralizedStirlingNumber_HPP_
#define GeneralizedStirlingNumber_HPP_
#include<cmath>
#include<climits>
#include<cfloat>
#include<cassert>
#include<unordered_map>
#include"util.hpp"
namespace AliasLDA {
class NMPair {
public:
	int n;
	int m;
	NMPair() :
			n(0), m(0) {
	}
	NMPair(int _n, int _m) :
			n(_n), m(_m) {
	}
};
class GeneralizedStirlingNumber {
	GeneralizedStirlingNumber() :
			lock(0) {
	}
	GeneralizedStirlingNumber(GeneralizedStirlingNumber const&);
	void operator=(GeneralizedStirlingNumber const&);

public:
	static GeneralizedStirlingNumber& getInstance() {
		static GeneralizedStirlingNumber instance;
		return instance;
	}
	std::unordered_map<int, double*> cache;
	std::unordered_map<int, NMPair> limits;
	std::unordered_map<int, int> reading;
	std::unordered_map<int, bool> writing;
	bool lock;
	inline int nm(int n, int m, int mLimit) {
		return n * mLimit + m;
	}
	double* getCache(double _a) {
		while (lock) {
		}
		lock = true;
		double* ret = NULL;
		int a = _a * 100;
		if (writing.count(a) == 0) {
			reading[a] = 1;
			ret = cache[a];
		} else {
			if (!writing[a]) {
				reading[a]++;
				ret = cache[a];
			}
		}
		lock = false;
		return ret;
	}
	void releaseCache(double _a) {
		while (lock) {
		}
		lock = true;
		int a = _a * 100;
		reading[a]--;
		lock = false;
	}
	NMPair getLimit(double _a) {
		while (lock) {
		}
		lock = true;
		NMPair ret;
		int a = _a * 100;
		if (writing.count(a) == 0) {
			//reading[a]=true;
			ret = limits[a];
		} else {
			if (!writing[a]) {
				//reading[a]=true;
				ret = limits[a];
			} else {
				ret.m = -1;
				ret.n = -1;
			}
		}
		lock = false;
		return ret;
	}
	bool initialize(double __a, int n, int m) {
		while (lock) {
		}
		lock = true;
		int a = __a * 100;
		if (reading.count(a) > 0) {
			if (reading[a]) {
				lock = false;
				return false;
			}
		}
		if (writing.count(a) > 0) {
			if (writing[a]) {
				lock = false;
				return false;
			}
		}
		writing[a] = true;
		if (cache.count(a) > 0) {
			NMPair pair = limits[a];
			int newn = std::max(n, pair.n);
			int newm = std::max(m, pair.m);
			if (newn != n && newm != m) {
				freeRecord(a);
			} else {
				writing[a] = false;
				lock = false;
				return true;
			}
		}
		//printf("test1\n");
		lock = false;
		NMPair p;
		p.n = n;
		p.m = m;
		limits[a] = p;
		double* localCache = new double[nm(n, 0, m)];
		//printf("test2\n");
		cache[a] = localCache;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				int ind = nm(i, j, m);
				localCache[ind] = -DBL_MAX;
			}
		}
		//printf("test3\n");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				int ind = nm(i, j, m);
				localCache[ind] = computelog(__a, i, j, m, localCache);
				//if(!isnan(localCache[ind]))printf("%f %d %d %d\n",localCache[ind],i,j,m);
			}
		}
		//printf("test4\n");
		writing[a] = false;
		return true;
	}
	~GeneralizedStirlingNumber() {
		//freeAll();
	}
	inline void freeRecord(int _a) {
		delete[] cache[_a];
		cache.erase(_a);
		limits.erase(_a);
	}
	inline void freeAll() {
		for (std::unordered_map<int, double*>::iterator it = cache.begin();
				it != cache.end(); it++) {
			delete[] it->second;
		}
	}
	inline double fact(double a, int n) {
		assert(n >= 0);
		if (n == 1) {
			return 0;
		}
		return (double) log((n - 1) - a) + fact(a, n - 1);
	}
	inline double computelog(double a, int n, int m, int mLimit,
			double* localCache) {
		if (m > n)
			return NAN;
		else if (m == n)
			return 0;
		else if (m <= 0)
			return n == 0 ? 0 : NAN;
		else {
			int ind = nm(n, m, mLimit);
			if (localCache[ind] != -DBL_MAX) {
				return localCache[ind];
			} else {
				double r;
				if (m == 1) {
					r = fact(a, n);
				} else {
					r = log(
							exp(
									computelog(a, n - 1, m - 1, mLimit,
											localCache)
											- computelog(a, n - 1, m, mLimit,
													localCache)) + (n - 1)
									- m * a)
							+ computelog(a, n - 1, m, mLimit, localCache);
				}
				localCache[ind] = r;
				return r;
			}
		}
	}

};
}
#endif
