///*
// * AliasTable.cpp
// *
// *  Created on: 15/01/2014
// *      Author: polymorpher
// */
//#include "util.hpp"
//#include"AliasTable.hpp"
//#include<sys/time.h>
//
//int main(){
//	VecInt labels;
//	labels.push_back(10);
//	labels.push_back(9);
//	labels.push_back(8);
//	labels.push_back(7);
//	labels.push_back(6);
//	labels.push_back(5);
//	labels.push_back(4);
//	labels.push_back(3);
//	labels.push_back(2);
//	labels.push_back(1);
//	AliasLDA::AliasTable atable(10,labels);
//	VecDouble probs;
//	probs.push_back(0.01);
//	probs.push_back(0.19);
//	probs.push_back(0.06);
//	probs.push_back(0.06);
//	probs.push_back(0.06);
//	probs.push_back(0.06);
//	probs.push_back(0.06);
//	probs.push_back(0.3);
//	probs.push_back(0.15);
//	probs.push_back(0.05);
//	atable.build(probs);
//
//	google::dense_hash_map<int,int> samples;
//	samples.set_empty_key(-1);
//			struct timeval tv;
//			struct timeval start_tv;
//			double elapsed = 0.0;
//			gettimeofday(&start_tv, NULL);
//	double iters=10000;
//	for(int i=0;i<iters;i++){
//		samples[atable.sample()]++;
//	}
//			gettimeofday(&tv, NULL);
//			elapsed = (tv.tv_sec - start_tv.tv_sec)
//					+ (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
//	for(auto it=samples.begin();it!=samples.end();it++){
//		printf("%d %lf\n",it->first,it->second / iters);
//	}
//	printf("atable size 10 sample with label test timespent=%lf avg=%lf\n",elapsed,elapsed/iters);
//
//	VecDouble prob2;
//	prob2.resize(1024,1/(double)1024);
//	AliasLDA::AliasTable atable2(1024);
//	atable2.build(prob2);
//	int csum2=0;
//	gettimeofday(&start_tv, NULL);
//	for(int i=0;i<iters;i++){
//		csum2+=atable2.sample();
//	}
//			gettimeofday(&tv, NULL);
//			elapsed = (tv.tv_sec - start_tv.tv_sec)
//					+ (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
//
//	printf("atable size 1024 sample test timespent=%lf avg=%lf csum2 avg=%lf\n",elapsed,elapsed/iters, csum2/iters);
//
//
//
//	VecDouble prob3;
//	prob3.resize(1024,1/(double)1024);
//	std::vector<AliasLDA::AliasTable> atables;
//	AliasLDA::AliasTable atable3(1024);
//	atables.resize(iters,atable3);
//	gettimeofday(&start_tv, NULL);
//	for(int i=0;i<iters;i++){
//		atables[i].build(prob3);
//	}
//			gettimeofday(&tv, NULL);
//			elapsed = (tv.tv_sec - start_tv.tv_sec)
//					+ (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
//
//	printf("build time test timespent=%lf avg=%lf csum2 avg=%lf\n",elapsed,elapsed/iters);
//
//
//	std::mt19937_64 rgen;
//	std::uniform_real_distribution<double> u01(0,1);
//	double csum=0;
//	gettimeofday(&start_tv, NULL);
//	for(int i=0;i<iters;i++){
//		csum+=u01(rgen);
//	}
//			gettimeofday(&tv, NULL);
//			elapsed = (tv.tv_sec - start_tv.tv_sec)
//					+ (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
//
//	printf("mt19937_64 test timespent=%lf avg=%lf d=%lf\n",elapsed,elapsed/iters,csum);
//	return 0;
//}
//
