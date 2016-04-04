/*
 * AliasLDA.hpp
 *
 *  Created on: 25/12/2013
 *      Author: polymorpher
 */

#pragma once

#include "ILDA.hpp"
#include <mutex>
#include "AliasTable.hpp"
#include<glog/logging.h>
#include<boost/format.hpp>
#include<boost/algorithm/string.hpp>
using boost::format;
namespace AliasLDA {

class AliasLDA_UltraSpeed: public AliasLDA::ILDA {
	friend class LDA;
	//priors
	double beta;
	double alpha;
	//rand
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;
	//Non-uniform Prior
	//VecDouble betaVec;
	VecDouble alphaVec;
	//Common between unifrom and non-unifrom prior
	double betasum;
	double alphasum;
	//basic
	VecInt nk;
	VecInt nw;
	Vector2DInt docs;
	Vector2DInt z;	//Topic assignment
	Vector2DInt nwk;
	int K;	//numTopics
	int V;	//sizeVocabulary
	//results
	Vector2DDouble phi;
	Vector2DDouble theta;
	VecDouble thetasum;

	//Sparse Stuff
	Vector2DInt ndkSparse;
	Vector2DInt ndkIndex;
	//bitshifting
	int oneWithTopicBitShift;
	int topicBitShift;
	int topicMask;
	int auxillaryMask;

	//AliasStuff
	//temp buffer
	int sampleFlagCount[3];
	//bucket stuff
	VecDouble bucket_sparse_values;
	VecInt bucket_sparse_indices;
	Vector2DInt aliasSamples;
	Vector2DDouble aliasSampleProps;
	VecDouble bucket_dense_sum;
	//threads
	std::thread* thread_table;
	bool terminated;
	VecInt aliasTableLocks;
	std::vector<AliasTable*> aliasTables;
	const int METRO_HASTING_STEPS = 2;
	const int ALIAS_SAMPLES_KEPT_RATIO = 16;
	const bool USE_MIXED_SUPPLY_CONSUME_MECH=true;

public:
	AliasLDA_UltraSpeed() :
			K(0), V(0), beta(0.1)/*,alpha(0.1)*/, betasum(0), alphasum(0), u01(
					0, 1), terminated(false) {
		rgen.seed(time(NULL)+rand());
	}
	virtual ~AliasLDA_UltraSpeed(){
		terminated=true;
		thread_table->join();
	}
	inline double rand01() {
		return u01(rgen);
	}
	inline int randInt(int limit) {
		return (int) (rand01() * limit);
	}
	inline int getnwk(int w, int k){
		return nwk[w][k];
	}
	inline int getndk(int d,int k){
		int index=ndkIndex[d][k];
		return index==0?0:decodeCount(ndkSparse[d][index-1]);
	}
	void initialiseBitShifts(){
		topicBitShift=0;
		topicMask=0;
		auxillaryMask=0;
		oneWithTopicBitShift=1;
		while(oneWithTopicBitShift<=K-1){
			topicMask+=oneWithTopicBitShift;
			oneWithTopicBitShift=oneWithTopicBitShift<<1;
			topicBitShift++;
		}
		auxillaryMask=~topicMask;
//		printf("topicBitShift=%x topicMask=%x, auxillaryMask=%x, oneWithTopicBitShift=%x\n ",
//				topicBitShift, topicMask, auxillaryMask, oneWithTopicBitShift);
	}
	void initialiseWordCount(){
		for(int d=0;d<docs.size();d++){
			for(int l=0;l<docs[d].size();l++){
				nw[docs[d][l]]++;
			}
		}
	}
	void initialiseWordTopicCounter(){
		nwk.resize(V);
		for(int w=0;w<V;w++){
			nwk[w].resize(K);
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
	void initialise_bucket_dense_sum(){
		bucket_dense_sum.resize(V);
		for(int w=0;w<V;w++){
			bucket_dense_sum[w]=0;
			for(int k=0;k<K;k++){
				int val=getnwk(w,k);
				if(val==0){
					bucket_dense_sum[w]+=(alpha*beta)/(nk[k]+betasum);
				}else{
					bucket_dense_sum[w]+=alpha*(beta+val)/(nk[k]+betasum);
				}
			}
		}
	}
	void initialise_bucket_sparse(){
		bucket_sparse_values.resize(K);
		bucket_sparse_indices.resize(K);
	}
	void initialise_sampleFlagCount(){
		sampleFlagCount[0]=0;
		sampleFlagCount[1]=0;
		sampleFlagCount[2]=0;
	}
	void initialiseAliasTables() {
		aliasTables.resize(V);
		aliasSamples.resize(V);
		aliasSampleProps.resize(V);
		aliasTableLocks.resize(V);
		for (int w = 0; w < aliasSamples.size(); w++) {
			if (USE_MIXED_SUPPLY_CONSUME_MECH
					&& (nw[w] * ALIAS_SAMPLES_KEPT_RATIO
							> sizeof(AliasTable)
									+ sizeof(AliasTable::AliasElement) * K)) {
				aliasSamples[w].resize(1);
				aliasSampleProps[w].resize(1);
				aliasSamples[w][0] = -1;
				aliasSampleProps[w][0] = -1;
			} else {
				aliasSamples[w].resize(nw[w] * ALIAS_SAMPLES_KEPT_RATIO + 1);
				aliasSampleProps[w].resize(
						nw[w] * ALIAS_SAMPLES_KEPT_RATIO + 1);
			}
		}
		for (int i = 0; i < V; i++)
			computeAliasTable(i);
		thread_table = new std::thread(
				&AliasLDA::AliasLDA_UltraSpeed::thread_computeAliasTable, this);
	}
	void initialise() {
		initialiseBeta();
		initialiseAlpha();
		initialiseWordTopicCounter();
		initialiseAssignmentOnly();
	}
	virtual void initialiseAssignmentOnly(){
		thetasum.resize(K);
		nw.resize(V);
		initialiseWordCount();
		initialiseDocumentTopicCounter();
		initialiseBitShifts();
		initialise_bucket_sparse();
		initialise_sampleFlagCount();
		initialiseTopicsAndCounts();
		initialiseAliasTables();
		LOG(INFO)<<"initialised.";
	}
	void initialiseTopicsAndCounts() {
		z.resize(docs.size());
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
		initialise_bucket_dense_sum();
	}
	void initialiseBeta() {
		betasum = beta * V;
	}
	void initialiseAlpha() {
		alphasum = 0;
		if (alpha != 0) {
			alphaVec.resize(K);
			for (size_t i = 0; i < alphaVec.size(); i++) {
				alphaVec[i] = alpha;
			}
		}
		for (auto it = alphaVec.begin(); it != alphaVec.end(); it++) {
			alphasum += *it;
		}
	}
	inline void inc_ndk(int d, int k) {
		int index=ndkIndex[d][k]-1;
		if(index<0){
			ndkSparse[d].push_back(oneWithTopicBitShift+k);
			index=ndkSparse[d].size()-1;
			ndkIndex[d][k]=index+1;
		}else{
			incEncodedTopic(ndkSparse[d][index]);
		}
		oneStepSortInc(ndkIndex[d],ndkSparse[d],index);
	}
	inline bool dec_ndk(int d, int k) {
		int index=ndkIndex[d][k]-1;
		assert(index>=0);
		assert(ndkSparse[d][index]>=0);
		decEncodedTopic(ndkSparse[d][index]);
		assert(ndkSparse[d].size()>index);
		return oneStepSortDec(ndkIndex[d],ndkSparse[d],index);
	}
	inline void inc_nwk(int w, int k) {
		nwk[w][k]++;
		nk[k]++;
	}
	inline void dec_nwk(int w, int k) {
		nwk[w][k]--;
		nk[k]--;
		assert(nwk[w][k]>=0);
	}
	inline int encodeTopic(int topic,int count){
		return (count<<topicBitShift)+topic;
	}
	inline int decodeTopic(int composite){
		return composite&topicMask;
	}
	inline int decodeCount(int composite){
		return (composite&auxillaryMask)>>topicBitShift;
	}
	inline void incEncodedTopic(int& composite){
		composite+=oneWithTopicBitShift;
	}
	inline void decEncodedTopic(int& composite){
		composite-=oneWithTopicBitShift;
	}
	void printCompositeVec(VecInt& vec){
		for(int i=0;i<vec.size();i++){
			printf("[%d,%d] ",decodeTopic(vec[i]),decodeCount(vec[i]));
		}
		printf("\n");
	}
	void printIndexMap(MapIntInt& m){
		for(auto it=m.begin();it!=m.end();it++){
			printf("[%d => %d] ",it->first,it->second);
		}
		printf("\n");
	}
	void oneStepSortInc(VecInt& m, VecInt& vec, const int index) {

		int val = vec[index];
		int topic = decodeTopic(val);
		int count= decodeCount(val);
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
		m[swappedTopic] = index+1;
		m[topic] = iterIndex +1 ;
	}
	bool oneStepSortDec(VecInt& m,VecInt& vec, int index) {
		int val = vec[index];
		int topic = decodeTopic(val);
		int count = decodeCount(val);
		if (count == 0) {
//			printf("erasing index=%d topic=%d\n",index,topic);
//			printf("currContentIndex=");printIndexMap(m);
//			printf("currContentVec=");printCompositeVec(vec);
			int iterIndex=vec.size()-1;
			int newval=vec[iterIndex];
			vec[index]=newval;
			m[decodeTopic(newval)]=index+1;
			m[topic]=0;
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
			m[swappedTopic] = index+1;
			m[topic] = iterIndex+1;
			return false;
		}
	}
	void update_bucket_dense_sum(int w,int k, bool inc){
		if(inc){
			bucket_dense_sum[w]-=alpha*(beta+getnwk(w,k)-1)/(nk[k]-1+betasum);
			bucket_dense_sum[w]+=alpha*computelwk(w,k);
		}else{
			bucket_dense_sum[w]-=alpha*(beta+getnwk(w,k)+1)/(nk[k]+1+betasum);
			bucket_dense_sum[w]+=alpha*computelwk(w,k);
		}
	}
	void computeAliasTable(int w){
//
//		AliasTable* table=new AliasTable(K);
//		table->build(lwkcache0[w]);
//		AliasTable* oldTable=denseAliasTables[w];
//		denseAliasTables[w]=table;
//		if(!aliasTableLocks[w]==1){
//			delete oldTable;
//		}
		VecDouble cache;
		cache.resize(K);
		double bds = 0;
		for(int k=0;k<K;k++){
			double val=(beta+getnwk(w,k))/(nk[k]+betasum);
			cache[k]=val;
			bds += val *alpha;
		}
		bucket_dense_sum[w]= bds;

		AliasTable* table=new AliasTable(K);
		table->build(cache);
		//if(nw[w]*ALIAS_SAMPLES_KEPT_RATIO < sizeof(AliasTable)+sizeof(AliasTable::AliasElement)*K){

		VecInt& samples=aliasSamples[w];
		VecDouble& props=aliasSampleProps[w];

		if(samples[0]!=-1){
			for(int i=0;i<nw[w]*ALIAS_SAMPLES_KEPT_RATIO;i++){
				int sample=table->sample();
				samples[i+1]=sample;
				props[i+1]=table->getProportions()[sample];
			}
			samples[0]=0;
			delete table;
		}else{
//			std::shared_ptr<AliasTable> oldTable=aliasTables[w];
			while(aliasTableLocks[w]==1){
				LOG(INFO)<<"WARNING [Table building thread]: Contention while replacing table!";
			}
			aliasTableLocks[w]=1;
//			while(!aliasTableLocks[w]->try_lock()){printf("WARNING [Table building thread]: Contention while replacing table!\n");}
			AliasTable* oldTable=aliasTables[w];
			aliasTables[w]=table;
//			aliasTableLocks[w]->unlock();
			aliasTableLocks[w]=0;
			if(oldTable!=NULL){
				delete oldTable;
			}
		}
	}
	void thread_computeAliasTable(){
		while(!terminated){
//			printf("rebuilding aliastable...\n");
			for(int w=0;w<V;w++){
				computeAliasTable(w);
			}
//			printf("rebuilding aliastable complete.\n");
		}
	}

	inline double computelwk(int w,int k){
		return (beta+getnwk(w,k))/(nk[k]+betasum);
	}


//	inline void discreteSample(VecInt& indices, VecDouble& prop, int size, VecDouble& rolls, VecInt& samples, int sampleSize){
//		int ptr=0;
//		double sum=0;
//		for(int i=0;i<size;i++){
//			sum+=prop[i];
//			if(ptr>=sampleSize)return;
//			while(ptr<sampleSize && rolls[ptr]<=sum){
//				samples[ptr]=indices[i];
//				ptr++;
//			}
//		}
//	}
	int discreteSample(VecInt& indices, VecDouble& prop, int size, double roll){
		for(int i=0;i<size;i++){
			if(roll<=prop[i])return indices[i];
			else roll-=prop[i];
		}
		return -1;
	}
	void metroHastingSample(int d,int w,int k,int l, double _bucket_sparse_sum,
			int discreteSampleResult, int& sample, int& sampleType){
		//AliasTable* at=denseAliasTables[w];
		VecInt& samples=aliasSamples[w];
		VecDouble& sampleProps=aliasSampleProps[w];
		int index=samples[0]+1;
		int ret=k;
//		double zprop=zSampleProp[d][l];
		//zprop= (zprop==0?computelwk(w,k):zprop);

		double sump=_bucket_sparse_sum+bucket_dense_sum[w];
		bool hasAccept=false;
		while(aliasTableLocks[w]==1){
			LOG(INFO)<<("WARNING: Contention inside core loop");
		}
		aliasTableLocks[w]=1;
//		while(!aliasTableLocks[w]->try_lock()){printf("WARNING: Contention inside core loop\n");};
		for(int i=0;i<METRO_HASTING_STEPS;i++){
			double roll0=rand01()*sump;
			bool accepted=false;
			if(roll0<=_bucket_sparse_sum){
				sample=discreteSampleResult;
				assert(sample!=-1);
				sampleType=0;
				accepted=true;
//				printf("*MH=%d d=%d l=%d w=%d k=%d sample=%d "
//						"type=%d sparse_sum=%lf roll=%lf\n",i,d,l,w,k,sample,
//						sampleType,bucket_sparse_sum,roll0);
			}else{
				double newZSampleProp;
				if(index!=0){
					sample=samples[index];
					newZSampleProp=sampleProps[index];
					index++;
				}else{
					AliasTable* table=aliasTables[w];
					assert(table!=NULL);
					sample=table->sample();
					newZSampleProp=table->getProportions()[sample];
				}
				sampleType=1;
				double propPrevSample=computelwk(w,ret);
				double propNewSample=computelwk(w,sample);
				int ndk0=getndk(d,ret);
				int ndk1=getndk(d,sample);
				double acceptance=((ndk1+alpha)/(ndk0+alpha))*(propNewSample/propPrevSample)*
						(( alpha*propPrevSample+propPrevSample*ndk0)/
						(alpha*newZSampleProp+propNewSample*ndk1));
				double roll=rand01();
				if(roll<acceptance){
					accepted=true;
				}
//				printf("MH=%d d=%d l=%d w=%d oldsample=%d sample=%d "
//						"type=%d sparse_sum=%lf ndk0=%d "
//						"ndk1=%d newZSampleProp=%lf propPrevSample=%lf propNewSample=%lf acceptance=%.13lf roll=%lf accepted=%d\n",
//						i,d,l,w,ret,sample,
//						sampleType, bucket_sparse_sum, ndk0,
//						ndk1,newZSampleProp,propPrevSample,propNewSample,acceptance,roll,accepted);

			}
			if(accepted){
				ret=sample;
				hasAccept=true;
			}
			assert(index<=nw[w]*ALIAS_SAMPLES_KEPT_RATIO);
		}
		aliasTableLocks[w]=0;
//		aliasTableLocks[w]->unlock();
		if(!hasAccept)sampleType=2;
		if(index!=0)samples[0]=index-1;
		sample=ret;
	}
	virtual void sampleWord(int d, int l) {
		int w = docs[d][l];
		int k = z[d][l];
		dec_ndk(d, k);
		dec_nwk(w, k);
		update_bucket_dense_sum(w,k,false);
		double bucket_sparse_sum=0;
		VecInt& ndkSparseLocal=ndkSparse[d];
		for(int i=0;i<ndkSparseLocal.size();i++){
			int encodedValue=ndkSparseLocal[i];
			int topic=decodeTopic(encodedValue);
			int count=decodeCount(encodedValue);
			double val=count*computelwk(w,topic);
			bucket_sparse_values[i]=val;
			bucket_sparse_indices[i]=topic;
			bucket_sparse_sum+=val;
		}
		int discreteSampleResult=discreteSample(bucket_sparse_indices,bucket_sparse_values,ndkSparseLocal.size(),rand01()*bucket_sparse_sum);
//		sparsePartTable->build(bucket_sparse_values);
		int newTopic=-1,flag=-1;
		metroHastingSample(d,w,k,l,bucket_sparse_sum,discreteSampleResult,newTopic, flag);
		sampleFlagCount[flag]++;
//		if(newTopic==-1 || flag==-1){
//			printf(" d=%d,l=%d | flag=%d oldtopic=%d newtopic=%d\n",d,l,flag ,k, newTopic);
//		}
		assert(newTopic != -1);
		inc_ndk(d, newTopic);
		inc_nwk(w, newTopic);
		update_bucket_dense_sum(w,newTopic,true);
		z[d][l] = newTopic;
	}
	virtual void gibbsStep() {
		for (size_t d = 0; d < docs.size(); d++) {
			VecInt& doc = docs[d];
			for (size_t l = 0; l < doc.size(); l++) {
				sampleWord(d, l);
			}
		}
		for(int i=0;i<3;i++){
			LOG(INFO)<<format("flags %d=%d") % i % sampleFlagCount[i];
			sampleFlagCount[i]=0;
		}
	}
	virtual void computePhi() {
		phi.resize(V);
		for (int w = 0; w < V; w++) {
			phi[w].resize(K);
			for(int k=0;k<K;k++){
				phi[w][k]=(nwk[w][k]+beta)/(double)(nk[k]+betasum);
			}
		}
	}
	virtual void computeTheta() {
		theta.resize(docs.size());
		for (size_t d = 0; d < docs.size(); d++) {
			theta[d].resize(K);
			for(int k=0;k<K;k++){
				theta[d][k]=(getndk(d,k)+alphaVec[k])/(double(docs[d].size()+alphasum));;

			}
		}
	}
	virtual void computeThetaSum() {
		for (int k = 0; k < K; k++) {
			thetasum[k] = 0;
		}
		for (size_t d = 0; d < docs.size(); d++) {
			VecInt flags;
			flags.resize(K);
			for (auto it = ndkSparse[d].begin(); it != ndkSparse[d].end();
					it++) {
				int encodedValue = *it;
				int topic = decodeTopic(encodedValue);
				int count = decodeCount(encodedValue);
				flags[topic] = 1;
				double thetaval = (count + alphaVec[topic])
						/ (double(docs[d].size() + alphasum));
				thetasum[topic] += thetaval / (double) docs.size();
			}
			for (int i = 0; i < K; i++) {
				if (flags[i] != 1) {
					double val = alphaVec[i] / (double) (docs[d].size() + alphasum);
					thetasum[i] += val / (double) docs.size();
				}
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
		for(int i=0;i<phi.size();i++){
			phi[i].clear();
			phi[i].shrink_to_fit();
		}
		phi.clear();
		phi.shrink_to_fit();
	}
	virtual void clearTheta(){
		for(int i=0;i<theta.size();i++){
			theta[i].clear();
			theta[i].shrink_to_fit();
		}
		theta.clear();
		theta.shrink_to_fit();
	}
	virtual void copyState(ILDA* _lda){
		AliasLDA_UltraSpeed* lda = static_cast<AliasLDA_UltraSpeed* > ( _lda);
		K = lda->K;
		V = lda->V;
		nk = lda->nk;
		nwk = lda->nwk;
		beta = lda->beta;
		alpha = lda->alpha;
		alphasum=lda->alphasum;
		alphaVec = lda->alphaVec;
		betasum = lda->betasum;
	}
	virtual void copyState(void* _lda, std::string impl){
		AliasLDA_UltraSpeed* lda = static_cast<AliasLDA_UltraSpeed* > ( _lda);
		K = lda->K;
		V = lda->V;
		nk = lda->nk;
		nwk = lda->nwk;
		beta = lda->beta;
		alpha = lda->alpha;
		alphasum=lda->alphasum;
		alphaVec = lda->alphaVec;
		betasum = lda->betasum;
	}

};

}
