/*
 * SparseLDA.hpp
 *
 *  Created on: 25/12/2013
 *      Author: polymorpher
 */

#pragma once

#include "ILDA.hpp"
#include<glog/logging.h>
#include<boost/format.hpp>
#include<boost/algorithm/string.hpp>
using boost::format;
namespace AliasLDA {

class SparseLDA_Speed: public AliasLDA::ILDA {
	friend class LDA;
	//basic
	VecInt nk;
	VecInt nw;
	Vector2DInt docs;
	Vector2DInt z;	//Topic assignment

	int K;	//numTopics
	int V;	//sizeVocabulary
	//bitshifting
	int oneWithTopicBitShift;
	int topicBitShift;
	int topicMask;
	int auxillaryMask;
	Vector2DInt nwkSparse;//Awesome SparseLDA special data structure
	Vector2DInt nwkIndex; //Indexing nwkSparse
	VecInt nwkLimits;
	Vector2DInt ndkSparse;
	Vector2DInt ndkIndex;
	VecInt ndkLimits;
	//bucket stuff
	double bucket_s;
	double bucket_r;
	double bucket_q;
	VecDouble bucket_s_values;
	VecDouble bucket_r_values;
	VecInt bucket_r_indices;
	int bucket_r_index_size;
	VecInt bucket_r_indices_indices;

	VecDouble bucket_q_values;
	VecDouble bucket_q_coeff;
	VecDouble bucket_q_coeff_enhanced;
	VecInt bucket_q_indices;
	int bucket_q_index_size;



	//priors
	double beta;
	double alpha;
	//Non-uniform Prior
	//VecDouble betaVec;
	VecDouble alphaVec;
	//Common between unifrom and non-unifrom prior
	double betasum;
	double alphasum;

	//rand
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;

	//results
	Vector2DDouble phi;
	Vector2DDouble theta;
	VecDouble thetasum;


	//temp buffer
	int sampleFlagCount[3];

public:
	SparseLDA_Speed() :
			K(0), V(0), beta(0.1)/*,alpha(0.1)*/, betasum(0), alphasum(0), u01(
					0, 1) {
		rgen.seed(time(NULL)+rand());
	}
	inline double rand01() {
		return u01(rgen);
	}
	inline int randInt(int limit) {
		return (int) (rand01() * limit);
	}
	void initialiseWordCount(){
		for(int d=0;d<docs.size();d++){
			for(int l=0;l<docs[d].size();l++){
				nw[docs[d][l]]++;
			}
		}
	}
	void initialiseWordTopicCounter(){
		nwkSparse.resize(V);
		nwkIndex.resize(V);
		for(int w=0;w<V;w++){
//			nwkSparse[w].resize(std::min(V,nw[w]));
			nwkSparse[w].reserve(std::min(V,nw[w]));
			nwkIndex[w].resize(K);
			for(int k=0;k<K;k++){
				nwkIndex[w][k]=0;
			}
		}
	}

	void initialiseDocumentTopicCounter(){
		ndkSparse.resize(docs.size());
		ndkIndex.resize(docs.size());
		for(int d=0;d<docs.size();d++){
			//ndkSparse[d].resize(std::min((int)docs[d].size(),K));
			ndkSparse[d].reserve(std::min((int)docs[d].size(),K));
			ndkIndex[d].resize(K);
			for(int k=0;k<K;k++){
				ndkIndex[d][k]=0;
			}
		}
	}
	void initialise() {
		initialiseBeta();
		initialiseAlpha();
		nw.resize(V);
		initialiseWordCount();
		initialiseWordTopicCounter();
		initialiseAssignmentOnly();
	}
	virtual void initialiseAssignmentOnly(){
		initialiseBitShifts();
		initialiseDocumentTopicCounter();
		thetasum.resize(K);
		bucket_q_coeff.resize(K);
		bucket_q_coeff_enhanced.resize(K);
		bucket_q_values.resize(K);
		bucket_s_values.resize(K);
		bucket_q_indices.resize(K);
		bucket_r_values.resize(K);
		bucket_r_indices.resize(K);
		bucket_r_indices_indices.resize(K);
		sampleFlagCount[0]=0;
		sampleFlagCount[1]=0;
		sampleFlagCount[2]=0;
		initialiseTopicsAndCounts();
	}
	void initialiseTopicsAndCounts() {
		z.resize(docs.size());
		LOG(INFO)<<"K="<<K;
		nk.resize(K);
		for (size_t d = 0; d < docs.size(); d++) {
			VecInt& doc = docs[d];
			z[d].resize(doc.size());
			VecInt& zd = z[d];
			for (size_t l = 0; l < doc.size(); l++) {
				int w = doc[l];
				int k = randInt(K);
				inc_ndk(d, k);
				inc_nwk(w, k);
				zd[l] = k;
			}
		}
	}
	void initialiseBeta() {
		betasum = beta * V;
	}
	void initialiseAlpha() {
		alphasum = alpha*K;
		alphaVec.resize(K);
		for(int i=0;i<alphaVec.size();i++){
			alphaVec[i]=alpha;
		}

	}
	inline int getnwk(int w, int k){
		int index=nwkIndex[w][k];
		return index<=0?0:decodeCount(nwkSparse[w][index-1]);
	}
	inline int getndk(int d,int k){
		int index=ndkIndex[d][k];
		return index<=0?0:decodeCount(ndkSparse[d][index-1]);
	}
	void initialiseBitShifts() {
		topicBitShift = 0;
		topicMask = 0;
		auxillaryMask = 0;
		oneWithTopicBitShift = 1;
		while (oneWithTopicBitShift <= K - 1) {
			topicMask += oneWithTopicBitShift;
			oneWithTopicBitShift = oneWithTopicBitShift << 1;
			topicBitShift++;
		}
		auxillaryMask = ~topicMask;
		//		printf("topicBitShift=%x topicMask=%x, auxillaryMask=%x, oneWithTopicBitShift=%x\n ",
		//				topicBitShift, topicMask, auxillaryMask, oneWithTopicBitShift);
	}
	inline void inc_ndk(int d, int k) {
		int index = ndkIndex[d][k] - 1;
		if (index < 0) {
			ndkSparse[d].push_back(oneWithTopicBitShift | k);
			index = ndkSparse[d].size() - 1;
			ndkIndex[d][k] = index + 1;
		} else {
			incEncodedTopic (ndkSparse[d][index]);
		}
		oneStepSortInc(ndkIndex[d], ndkSparse[d], index);
	}
	inline bool dec_ndk(int d, int k) {
		int index = ndkIndex[d][k] - 1;
		CHECK(index >= 0);
		CHECK(ndkSparse[d][index] >= 0);
		decEncodedTopic (ndkSparse[d][index]);
		CHECK(ndkSparse[d].size() > index);
		return oneStepSortDec(ndkIndex[d], ndkSparse[d], index);
	}
	inline void inc_nwk(int w, int k) {
		int index = nwkIndex[w][k] - 1;
		if (index < 0) {
			nwkSparse[w].push_back(oneWithTopicBitShift | k);
			index = nwkSparse[w].size() - 1;
			nwkIndex[w][k] = index + 1;
		} else {
			incEncodedTopic (nwkSparse[w][index]);
		}
		nk[k]++;
		oneStepSortInc(nwkIndex[w], nwkSparse[w], index);

	}
	inline bool dec_nwk(int w, int k) {
		int index = nwkIndex[w][k] - 1;
		CHECK(index != -1);
		CHECK(nwkSparse[w][index] >= 0);
		decEncodedTopic (nwkSparse[w][index]);
		nk[k]--;
		return oneStepSortDec(nwkIndex[w], nwkSparse[w], index);

	}
	inline int encodeTopic(int topic, int count) {
		return (count << topicBitShift) | topic;
	}
	inline int decodeTopic(int composite) {
		return composite & topicMask;
	}
	inline int decodeCount(int composite) {
		return (composite & auxillaryMask) >> topicBitShift;
	}
	inline void incEncodedTopic(int& composite) {
		composite += oneWithTopicBitShift;
	}
	inline void decEncodedTopic(int& composite) {
		composite -= oneWithTopicBitShift;
	}
	void printCompositeVec(VecInt& vec) {
		for (int i = 0; i < vec.size(); i++) {
			printf("[%d,%d] ", decodeTopic(vec[i]), decodeCount(vec[i]));
		}
		printf("\n");
	}

	void oneStepSortInc(VecInt& m, VecInt& vec, const int index) {

		int val = vec[index];
		int topic = decodeTopic(val);
		int count = decodeCount(val);
		int iterIndex = index - 1;
		//printf("iterindex=%d\n",iterIndex);
		while (iterIndex >= 0) {
			if (decodeCount(vec[iterIndex]) >= count) {
				break;
			}
			iterIndex--;
		}
		iterIndex++;
		//printf("iterindex=%d\n",iterIndex);
		int swappedTopic = decodeTopic(vec[iterIndex]);
		vec[index] = vec[iterIndex];
		vec[iterIndex] = val;
		m[swappedTopic] = index + 1;
		m[topic] = iterIndex + 1;
	}
	bool oneStepSortDec(VecInt& m, VecInt& vec, int index) {
		CHECK(index < vec.size());
		int val = vec[index];
		int topic = decodeTopic(val);
		int count = decodeCount(val);
		if (count == 0) {
			//			printf("erasing index=%d topic=%d\n",index,topic);
			//			printf("currContentIndex=");printIndexMap(m);
			//			printf("currContentVec=");printCompositeVec(vec);
			int iterIndex = vec.size() - 1;
			int newval = vec[iterIndex];
			vec[index] = newval;
			m[decodeTopic(newval)] = index + 1;
			m[topic] = 0;
			vec[iterIndex]=0;
			vec.resize(iterIndex);
			//			printf("afterContentIndex=");printIndexMap(m);
			//			printf("afterContentVec=");printCompositeVec(vec);
			return true;
		} else {
			int iterIndex = index + 1;
			while (iterIndex < vec.size()) {
				if (decodeCount(vec[iterIndex]) <= count) {
					break;
				}
				iterIndex++;
			}
			iterIndex--;
			int swappedTopic = decodeTopic(vec[iterIndex]);
			vec[index] = vec[iterIndex];
			vec[iterIndex] = val;
			m[swappedTopic] = index + 1;
			m[topic] = iterIndex + 1;
			return false;
		}
	}

	void computeBucketS(){
		bucket_s=0;
		for(int k=0;k<K;k++){
			double val=(alpha*beta)/(double)(betasum+nk[k]);
			bucket_s+=val;
			bucket_s_values[k]=val;
		}
	}
	inline void updateBucketS(int k){
		bucket_s-=bucket_s_values[k];
		double newval=(alpha*beta)/(double)(betasum+nk[k]);
		bucket_s+=newval;
		bucket_s_values[k]=newval;
	}
	void computeBucketR(int d){
		bucket_r=0;
		bucket_r_index_size=0;
		for(int i=0;i<K;i++){
			bucket_r_indices_indices[i]=0;
		}
		//for(auto it=ndkSparse[d].begin();it!=ndkSparse[d].end();it++){
		VecInt& ndkSparseLocal=ndkSparse[d];
		for(int i=0;i<ndkSparseLocal.size();i++){
			//int encodedValue=*it;
			int encodedValue=ndkSparseLocal[i];
			int topic=decodeTopic(encodedValue);
			int count=decodeCount(encodedValue);
			CHECK(!(topic<=0&&count<=0));
			double val=(count*beta)/(double)(betasum+nk[topic]);
			bucket_r+=val;
			bucket_r_values[topic]=val;
			bucket_r_indices[bucket_r_index_size]=topic;
			bucket_r_indices_indices[topic]=bucket_r_index_size+1;
//			printf("%d bucket_r oval=%d topic=%d count=%d val=%lf bucket_r_index_size=%d\n"
//					,d,encodedValue,topic,count,val,bucket_r_index_size);
			bucket_r_index_size++;

		}
	}
	inline void updateBucketR(int d,int k, bool removed_ndk){
		int index = bucket_r_indices_indices[k] - 1;
//		CHECK(index>=0);
		if (index < 0 ) {
			double newval = (getndk(d,k) * beta) / (double) (betasum + nk[k]);
			CHECK(!removed_ndk);
			bucket_r += newval;
//			printf("update %d %d, 1add %lf\n", d, k, newval);
			CHECK(bucket_r_index_size < K);
			bucket_r_values[k]=newval;
			bucket_r_indices[bucket_r_index_size] = k;
			bucket_r_indices_indices[k] = bucket_r_index_size + 1;
			bucket_r_index_size++;
			CHECK(bucket_r_index_size<=K);
		} else {
			if(bucket_r_indices[index]!=k){
				LOG(INFO)<<format("%1% bucket_r_indices[index]=%2%, k=%3%") % d % bucket_r_indices[index] % k;
				CHECK(bucket_r_indices[index]==k);
			}
			bucket_r -= bucket_r_values[k];
//			printf("update %d %d, subtracting %lf\n", d, k, bucket_r_values[k]);
//			CHECK(bucket_r >= 0.0);
			if(bucket_r<=0){
//				LOG(INFO)<<"bucket_r=="<<bucket_r;
				bucket_r=0;
			}
			if (!removed_ndk) {
				double newval = (getndk(d,k) * beta) / (double) (betasum + nk[k]);
				bucket_r += newval;
				bucket_r_values[k] = newval;
//				printf("update %d %d, 2add %lf count=%d\n", d, k, newval,count);
			} else {
//				printf("update %d %d, removing %lf\n", d, k, bucket_r_values[k]);
				bucket_r_values[k] = 0;
				double newTopic=bucket_r_indices[bucket_r_index_size - 1];
				bucket_r_indices[index]=newTopic;
				bucket_r_indices[bucket_r_index_size - 1]=-1;
				bucket_r_indices_indices[newTopic] = index+1;
				bucket_r_indices_indices[k] = 0;
				CHECK(bucket_r_index_size>=1);
				bucket_r_index_size--;

			}
		}
	}
	void computeBucketQCoeffConst(){
		for(int k=0;k<K;k++){
			double val=alpha/(double)(betasum+nk[k]);
			bucket_q_coeff[k]=val;
			bucket_q_coeff_enhanced[k]=val;
			//printf("computeBucketQCoeff: bucket_q_coeff_enhanced[%d]=%lf\n",k,val);
		}
	}
	void computeBucketQCoeffRemaining(int d){
		//for(auto it=ndkSparse[d].begin();it!=ndkSparse[d].end();it++){
		VecInt& ndkLocal=ndkSparse[d];
		for(int i=0;i<ndkLocal.size();i++){
			//int encodedValue=*it;
			int encodedValue=ndkLocal[i];
			int topic=decodeTopic(encodedValue);
			int count=decodeCount(encodedValue);
			CHECK(!(topic==0&&count==0));
			bucket_q_coeff_enhanced[topic]+=count/(double)(betasum+nk[topic]);
			//printf("ndksparse[%d][%d]=%d,bucket_q_coeff_enhanced[%d]=%lf\n",d,k,val,k,bucket_q_coeff_enhanced[k]);
		}
	}
	void resetBucketQCoeffRemaining(int d){
		for(int k=0;k<K;k++){
			bucket_q_coeff_enhanced[k]=bucket_q_coeff[k];
		}
	}
	inline void updateBucketQCoeffAll(int d, int k, bool removed_ndk){

		double val=alpha/(double)(betasum+nk[k]);
		bucket_q_coeff[k]=val;
		if(removed_ndk){
			bucket_q_coeff_enhanced[k]=val;
		}else{
			bucket_q_coeff_enhanced[k]=val + (getndk(d,k)/(double)(betasum+nk[k]));
		}
	}
	void computeBucketQValues(int d, int w){
		VecInt& nwkSparseLocal=nwkSparse[w];
		bucket_q_index_size=0;
		bucket_q=0;
		//for(auto it=nwkSparseLocal.begin();it!=nwkSparseLocal.end();it++){
		for(int i=0;i<nwkSparseLocal.size();i++){
			//int val=*it;
			int val=nwkSparseLocal[i];
			int topic=decodeTopic(val);
			int count=decodeCount(val);
			CHECK(!(topic==0&&count==0));
			double qval=bucket_q_coeff_enhanced[topic]*count;
//			printf("oval=%d, topic=%d, count=%d,bucket_q_coeff_enhanced[%d]=%lf, bucket_q_values[%d]=%lf\n"
//					,val,topic,count,topic,bucket_q_coeff_enhanced[topic],bucket_q_index_size,qval);
			bucket_q_indices[bucket_q_index_size]=topic;
			bucket_q_values[bucket_q_index_size]=qval;
			bucket_q_index_size++;
			bucket_q+=qval;
		}
	}



	int sampleFromBucketQ(int d, int l, double roll){
		int newTopic = -1;
		for (int i = 0; i < bucket_q_index_size; i++) {
			double prob = bucket_q_values[i];
//            double prob = 
			//				printf("prob=%lf\n",prob);
			if (roll <= prob) {
				newTopic = bucket_q_indices[i];
				break;
			} else {
				roll -= prob;
				//					printf("roll=%lf\n",roll);
			}
		}
		return newTopic;
	}
	int sampleFromBucketR(int d,int l,double roll){
		int newTopic=-1;
		//for(auto it=bucket_r_values.begin();it!=bucket_r_values.end();it++){
		for(int i=0;i<bucket_r_index_size;i++){
			int topic=bucket_r_indices[i];
			double prob=bucket_r_values[topic];
			if(roll<=prob){
				newTopic=topic;
				break;
			}else{
				roll-=prob;
			}
		}
		return newTopic;
	}
	int sampleFromBucketS(int d,int l,double roll){
		int newTopic=-1;
		for (int i = 0; i < K; i++) {
			double sval = bucket_s_values[i];
			if (roll <= sval) {
				newTopic = i;
				break;
			} else {
				roll -= sval;
			}
		}
		return newTopic;
	}
	virtual void sampleWord(int d, int l) {
		int w = docs[d][l];
		int k = z[d][l];
		bool removed_ndk=dec_ndk(d, k);
		bool removed_nwk=dec_nwk(w, k);
		updateBucketS(k);
		updateBucketR(d,k,removed_ndk);
		updateBucketQCoeffAll(d,k,removed_ndk);
		double sump = 0;
		if(removed_nwk){
//			VecInt& nwkSparseLocal=nwkSparse[w];
//			CHECK(getnwk(w,k)==0);
//			CHECK(nwkSparseLocal.size()==0);
			computeBucketQValues(d,w);
//			double vv=0.0;
//			if(bucket_q>0){
//				LOG(INFO)<<"bucket_q="<<bucket_q;
//				LOG(INFO)<<nwkSparseLocal.size();
//				CHECK(bucket_q<0.00000001L);
//			}

			sump=bucket_q+bucket_s+bucket_r;
		}else{
			computeBucketQValues(d,w);
			sump=bucket_s+bucket_r+bucket_q;
		}
		double roll = rand01() * sump;
		double originalRoll=roll;
		int newTopic = -1;
		int flag=-1;
//		if((!removed_nwk) &&roll<=bucket_q){
		if(roll<=bucket_q){
			flag=0;
			newTopic=sampleFromBucketQ(d,l,roll);
		}else{
			roll-=bucket_q;
			if(roll<=bucket_r){
				flag=1;
				newTopic=sampleFromBucketR(d,l,roll);
			}else{
				roll -= bucket_r;
				if (roll <= bucket_s) {
					flag = 2;
					newTopic=sampleFromBucketS(d,l,roll);
				}
			}
		}
		sampleFlagCount[flag]++;

//		if (abs(r) <= 0.0001 && newTopic == -1) {
//			newTopic = K - 1;
//		}
//		if(newTopic==-1 || flag==-1){
//			printf("removed_nwk=%d, d=%d,l=%d | flag=%d q=%.14lf r=%.14lf s=%.14lf roll=%.14lf oroll=%.14lf newtopic=%d\n",
//					removed_nwk,d,l,flag,bucket_q,bucket_r,bucket_s,roll, originalRoll,newTopic);
//		}
		CHECK(newTopic != -1);
		inc_ndk(d, newTopic);
		inc_nwk(w, newTopic);
		updateBucketS(newTopic);
		updateBucketR(d,newTopic,false);
		updateBucketQCoeffAll(d,newTopic,false);
		z[d][l] = newTopic;
	}
	virtual void gibbsStep() {
		computeBucketS();
		computeBucketQCoeffConst();
		for (size_t d = 0; d < docs.size(); d++) {
			computeBucketR(d);
			computeBucketQCoeffRemaining(d);
			VecInt& doc = docs[d];
			for (size_t l = 0; l < doc.size(); l++) {
				sampleWord(d, l);
			}
			resetBucketQCoeffRemaining(d);
		}
		for(int i=0;i<3;i++){
			LOG(INFO)<<format("flags %1%=%2%") % i % sampleFlagCount[i];
			sampleFlagCount[i]=0;
		}
	}
	virtual void computePhi() {
		phi.resize(V);

		for (int w = 0; w < V; w++) {
			phi[w].resize(K);
			for (int k = 0; k < K; k++) {
				phi[w][k] = (getnwk(w, k) + beta) / (double) (nk[k] + betasum);
			}
		}
	}
	virtual void computeThetaSum() {
		for (int k = 0; k < K; k++) {
			thetasum[k] = 0;
		}
		for (size_t d = 0; d < docs.size(); d++) {
			for (int k = 0; k < K; k++) {
				double val = (getndk(d, k) + alpha)
						/ (double) (docs[d].size() + alphasum);
				thetasum[k] += val / (double) docs.size();
			}
		}
	}
	virtual void computeTheta() {
		theta.resize(docs.size());
		for (size_t d = 0; d < docs.size(); d++) {
			theta[d].resize(K);
			for(int k=0;k<K;k++){
				theta[d][k]=(getndk(d,k)+alpha)/(double)(docs[d].size()+alphasum);

			}
		}
	}
	//getters and setters
	virtual void setDocuments(Vector2DInt& docs) {
		this->docs.resize(docs.size());
		for (size_t i = 0; i < docs.size(); i++) {
			auto& doc = this->docs[i];
			auto& thisdoc = docs[i];
			doc.resize(docs[i].size());
			for (size_t j = 0; j < thisdoc.size(); j++) {
				doc[j] = thisdoc[j];
				//printf("%d ",thisdoc[j]);
			}
		}
		//this->docs=docs;
	}
	//	double getAlpha() const {
	//		return alpha;
	//	}
	//
	//	void setAlpha(double alpha) {
	//		this->alpha = alpha;
	//	}

	virtual const VecDouble& getAlphaVec() const {
		return alphaVec;
	}

	virtual void setAlphaVec(const VecDouble& alphaVec) {
		this->alphaVec = alphaVec;
	}
	virtual void setAlpha(double alpha) {
		this->alpha = alpha;
	}
	virtual double getBeta() const {
		return beta;
	}

	virtual void setBeta(double beta) {
		this->beta = beta;
	}
	virtual void setNumTopics(int K) {
		this->K = K;
	}
	virtual void setSizeVocabulary(int V) {
		this->V = V;
	}
	virtual const Vector2DDouble& getPhi() const {
		return phi;
	}
	virtual const Vector2DDouble& getTheta() const {
		return theta;
	}
	virtual const VecDouble& getThetasum() const {
		return thetasum;
	}
	virtual int getNumTopics() const {
		return K;
	}
	virtual const Vector2DInt& getDocuments() const {
		return docs;
	}

	virtual void clearPhi(){
		phi.clear();
	}
	virtual void clearTheta(){
		theta.clear();
	}
	virtual void copyState(ILDA* _lda) {
		LOG(INFO)<<"Copying State...";
		SparseLDA_Speed* lda = static_cast<SparseLDA_Speed* > ( _lda);
		K = lda->K;
		V = lda->V;
		nk = lda->nk;
		nwkSparse = lda->nwkSparse;
		nwkIndex=lda->nwkIndex;
		beta = lda->beta;
		alpha = lda->alpha;
		alphasum=lda->alphasum;
		alphaVec = lda->alphaVec;
		betasum = lda->betasum;
	}
	virtual void copyState(void* _lda, std::string impl) {
		LOG(INFO)<<"Copying State...";
		SparseLDA_Speed* lda = static_cast<SparseLDA_Speed* > ( _lda);
		K = lda->K;
		V = lda->V;
		nk = lda->nk;
		nwkSparse = lda->nwkSparse;
		nwkIndex=lda->nwkIndex;
		beta = lda->beta;
		alpha = lda->alpha;
		alphasum=lda->alphasum;
		alphaVec = lda->alphaVec;
		betasum = lda->betasum;
	}
};

}
