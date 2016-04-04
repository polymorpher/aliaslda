/*
 * PerplexityAnalyzer.hpp
 *
 *  Created on: 04/01/2013
 *      Author: polymorpher
 */

#pragma once
#include"util.hpp"
inline bool isnan(double x){return x!=x;}
class PerplexityAnalyzer {
	int K;
	int V;
	const Vector2DInt& docs;
	const Vector2DDouble& phi;
	const Vector2DDouble& theta;
public:
	PerplexityAnalyzer(const Vector2DInt &_docs, const Vector2DDouble &_phi,
			const Vector2DDouble &_theta, int _K, int _V): docs(_docs),phi(_phi),theta(_theta),K(_K),V(_V){
	}
	double getPerplexity() {
		float sum = 0;
		int N = 0;
		for (int d = 0; d < docs.size(); d++) {
			N += docs[d].size();
			double dsum = 0;
			for (int l = 0; l < docs[d].size(); l++) {
				double wsum = 0;
				int w = docs[d][l];
				for(int k=0;k<K;k++){
					wsum+=theta[d][k]*phi[w][k];
					if(isnan(theta[d][k]) || isnan(phi[w][k])
							|| theta[d][k]<0 || phi[w][k]<0){
						printf("Warning  d=%d w=%d k=%d theta=%f phi=%f\n",d,w,k,theta[d][k],phi[w][k]);
					}
				}
				dsum += log(wsum);
			}
			sum += dsum;
		}
		printf("%lf %d\n", sum, N);
		return (double) exp(-sum / (double) N);
	}
};

