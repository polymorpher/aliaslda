
#pragma once

#include "GeneralizedStirlingNumber.hpp"
#include "ILDA.hpp"
#include"AliasTable.hpp"
namespace AliasLDA {
class AliasWordPDPLDA: public ILDA {
	friend class WordPDPLDA;
	const int METRO_HASTING_STEPS = 2;
	const int ALIAS_SAMPLES_KEPT_RATIO = 8;
	const bool USE_MIXED_SUPPLY_CONSUME_MECH=true;
	const double MIN_ERROR_BOUND=1e-7;

	Vector2DInt docs;
	int V;
	int K;
	double a; //discount
	double b; //concentration
	float gamma;  //for generating base Dirichlet
	float gammasum;
	VecDouble phi0;	//dim K V. base measure for each group. This variable is unused.
	Vector2DDouble phi; //dim K V. Word|topic,group distribution.
	Vector2DDouble theta; //dim I D K topic|document,group distribution.


	VecDouble thetasum; //dim I K theta_{I,D,K} sum over all documents
	double alpha;
	VecDouble alphaVec;
	double alphasum;
	Vector2DInt m; //dim V K. ==nwk;
	VecInt msum; //dim K.
	Vector2DInt t;  //dim V K. multiplicity
	VecInt tsum; //dim K
	Vector2DInt ndkSparse;
	Vector2DInt ndkIndex;
	std::vector<AliasTable*> aliasTables;
	Vector2DInt z; //dim I _D __L. topic assignment
	VecInt nw;
	int maxN;
	int maxM;
	int maxV;
	//bitshifting
	int oneWithTopicBitShift;
	int topicBitShift;
	int topicMask;
	int auxillaryMask;


	VecDouble bucket_sparse_values;
	VecInt bucket_sparse_indices;
	Vector2DInt aliasSamples;
	Vector2DDouble aliasSampleProps;
	VecDouble bucket_dense_sum;

	std::mt19937_64 rgen;
	std::uniform_real_distribution<float> u01;
	GeneralizedStirlingNumber& gsn;

	std::thread* thread_table;
	bool terminated;
	VecInt aliasTableLocks;
	int sampleFlagCount[3];
	int numMetroRej;

	Vector2DDouble gsnCache;
	int gsn_maxN;
	int gsn_maxM;
	void initialiseGSN() {
		NMPair limits = gsn.getLimit(a);
		while (limits.m == -1)
			limits = gsn.getLimit(a);
		gsn_maxN = limits.n;
		gsn_maxM = limits.m;
		if (limits.m < gsn_maxM || limits.n < gsn_maxN) {
			printf("initialising GSN Module...\n");
			while (!gsn.initialize(a, gsn_maxN, gsn_maxM)) {
			}
			printf("GSN Module initialized\n");
		}
		double* tempCache;
		gsnCache.resize(gsn_maxN);
		while ((tempCache = gsn.getCache(a)) == NULL) {
		}
		printf("Copying over GSN to local buffer\n");
		for (int i = 0; i < limits.n; i++) {
			gsnCache[i].resize(gsn_maxM);
			for (int j = 0; j < std::min(limits.m, i); j++) {
				int ind = i * limits.m + j;
				gsnCache[i][j] = tempCache[ind];
				//				printf("gsn %d %d val=%lf\n",i,j,gsnCache[i][j]);
			}
		}
	}
	double getGSNRatio(int n0, int m0, int n1, int m1) {
		assert(n0 < gsn_maxN && n1 < gsn_maxN);
		assert(m0 < gsn_maxM && m1 < gsn_maxM);
		assert(gsnCache[n0][m0]!=NAN && gsnCache[n1][m1]!=NAN);
		//		if(gsnCache[n0][m0] < gsnCache[n1][m1] ){
		//		printf("lS(%d,%d)=%lf, lS(%d,%d)=%lf, ratio=%lf\n",n0,m0,gsnCache[n0][m0],n1,m1,gsnCache[n1][m1],exp(gsnCache[n0][m0]-gsnCache[n1][m1]));
		//		}
		return exp(gsnCache[n0][m0] - gsnCache[n1][m1]);
	}

	inline void inc_m(int w, int k) {	//increment m
		m[w][k]++;
		msum[k]++;
	}
	inline void dec_m(int w, int k) {	//decrement m
		m[w][k]--;
		msum[k]--;
	}
	inline void inc_t(int w, int k) {
		t[w][k]++;
		tsum[k]++;
	}
	inline void dec_t(int w, int k) {
		t[w][k]--;
		tsum[k]--;
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
//		assert(ndkSparse[d].size()>index);
		return oneStepSortDec(ndkIndex[d],ndkSparse[d],index);
	}
	inline int encodeTopic(int topic, int count) {
		return (count << topicBitShift) + topic;
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
	void printIndexMap(MapIntInt& m) {
		for (auto it = m.begin(); it != m.end(); it++) {
			printf("[%d => %d] ", it->first, it->second);
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
		assert(vec.size()>index);
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
  	inline double getphi0(int w,int k){
  		return (gamma+t[w][k])/(gammasum+tsum[k]);
  	}


public:
  	AliasWordPDPLDA() :
			K(0), V(0), alphasum(0), u01(
					0, 1), terminated(false),gsn(GeneralizedStirlingNumber::getInstance()) {
		rgen.seed(time(NULL)+rand());
	}
	~AliasWordPDPLDA(){
		terminated=true;
		thread_table->join();
		delete thread_table;
	}
	inline double rand01() {
		return u01(rgen);
	}
	inline int randInt(int limit) {
		return (int) (rand01() * limit);
	}

	inline int getndk(int d, int k) {
		int index = ndkIndex[d][k];
		return index == 0 ? 0 : decodeCount(ndkSparse[d][index - 1]);
	}
	void initialiseBitShifts() {
		topicBitShift = 0;
		topicMask = 0;
		auxillaryMask = 0;
		oneWithTopicBitShift = 1;
		while (oneWithTopicBitShift < K - 1) {
			topicMask += oneWithTopicBitShift;
			oneWithTopicBitShift = oneWithTopicBitShift << 1;
			topicBitShift++;
		}
		auxillaryMask = ~topicMask;
	}
	void initialiseWordCount() {
		for (int d = 0; d < docs.size(); d++) {
			for (int l = 0; l < docs[d].size(); l++) {
				nw[docs[d][l]]++;
			}
		}
	}
	void initialiseRestaurantCounter() {
		m.resize(V);
		t.resize(V);
		msum.resize(K);
		tsum.resize(K);
		for (int w = 0; w < V; w++) {
			m[w].resize(K);
			t[w].resize(K);
		}
	}

	void initialiseDocumentTopicCounter() {
		ndkSparse.resize(docs.size());
		ndkIndex.resize(docs.size());
		for (int d = 0; d < docs.size(); d++) {
			//ndkSparse[d].resize(std::min((int)docs[d].size(),K));
			ndkSparse[d].reserve(std::min((int) docs[d].size(), K));
			ndkIndex[d].resize(K);
			for (int k = 0; k < K; k++) {
				ndkIndex[d][k] = 0;
			}
		}
	}
	inline double calculate_f_flat(int w, int composite) {
		return calculate_f(w, composite / 2, composite % 2 == 1);
	}
	inline double calculate_f_flat_withoffset(int w, int composite, int offset,
			int toffset) {
		return calculate_f_withoffset(w, composite / 2, composite % 2 == 1,
				offset, toffset);
	}
	inline double calculate_f(int w, int k, bool r) {
		if (r) {
			//			printf("calc f true\n");
			return ((b + a * tsum[k]) / (double) (b + msum[k]))
					* ((t[w][k] + 1) / (double) (m[w][k] + 1))
					* getGSNRatio(m[w][k] + 1, t[w][k] + 1, m[w][k], t[w][k])
					* ((gamma + t[w][k]) / (double) (gammasum + tsum[k]));
		} else {
			//			printf("calc f false\n");
			if (t[w][k] > 0 && m[w][k] > 0) {
				return ((m[w][k] - t[w][k] + 1) / (double) (m[w][k] + 1))
						* getGSNRatio(m[w][k] + 1, t[w][k], m[w][k], t[w][k])
						/ (double) (b + msum[k]);

			} else {
				return 0;
			}
		}
	}
	inline double calculate_f_withoffset(int w, int k, bool r, int offset,
			int toffset) {
		if (r) {
			return ((b + a * (tsum[k] + toffset))
					/ (double) (b + msum[k] + offset))
					* ((t[w][k] + 1 + toffset) / (double) (m[w][k] + 1 + offset))
					* getGSNRatio(m[w][k] + 1 + offset, t[w][k] + 1 + toffset,
							m[w][k] + offset, t[w][k] + toffset)
					* ((gamma + t[w][k] + toffset)
							/ (double) (gammasum + tsum[k] + toffset));
		} else {
			if (t[w][k] + toffset > 0 && m[w][k] + offset > 0) {
				return ((m[w][k] + offset - (t[w][k] + toffset) + 1)
						/ (double) (m[w][k] + offset + 1))
						* getGSNRatio(m[w][k] + offset + 1, t[w][k] + toffset,
								m[w][k] + offset, t[w][k] + toffset)
						/ (double) (b + msum[k] + offset);
			} else {
				return 0;
			}
		}
	}
	void initialise_bucket_dense_sum() {
		for (int w = 0; w < V; w++) {
			bucket_dense_sum[w] = 0;
			for (int k = 0; k < K; k++) {
				bucket_dense_sum[w] += alpha*calculate_f(w,k,true);
				bucket_dense_sum[w] += alpha*calculate_f(w,k,false);
			}
		}
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
	void initialiseTopicsAndCounts() {
		z.resize(docs.size());
		for (size_t d = 0; d < docs.size(); d++) {
			VecInt& doc = docs[d];
			z[d].resize(doc.size());
			VecInt& zd = z[d];
			for (size_t l = 0; l < doc.size(); l++) {
				int w = doc[l];
				int k = randInt(K);
				inc_ndk(d, k);
				inc_m(w, k);
				zd[l] = k;
				if(rand01()<=1/(double)(t[w][k]+1) || t[w][k]==0){
					inc_t(w,k);
				}
			}
		}
		initialise_bucket_dense_sum();
	}
	void initialiseGamma(){
		gammasum=gamma*V;
	}
	void initialise() {
		initialiseGSN();
		initialiseAlpha();
		initialiseGamma();
		initialiseRestaurantCounter();
		initialiseAssignmentOnly();
		printf("initialised.");
	}
	virtual void initialiseAssignmentOnly(){
		initialiseDocumentTopicCounter();
		initialiseBitShifts();
		nw.resize(V);
		initialiseWordCount();
		thetasum.resize(K);
		bucket_sparse_values.resize(2*K);
		bucket_sparse_indices.resize(2*K);
		bucket_dense_sum.resize(V);
		sampleFlagCount[0] = 0;
		sampleFlagCount[1] = 0;
		sampleFlagCount[2] = 0;
		initialiseTopicsAndCounts();
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
				&AliasLDA::AliasWordPDPLDA::thread_computeAliasTable, this);
	}
	void thread_computeAliasTable() {
		while (!terminated) {
			//			printf("rebuilding aliastable...\n");
			for (int w = 0; w < V; w++) {
				computeAliasTable(w);
			}
			//			printf("rebuilding aliastable complete.\n");
		}
	}
	void computeAliasTable(int w) {
		VecDouble cache;
		cache.resize(2*K);
		for (int k = 0; k < K; k++) {
			cache[2*k] = calculate_f(w, k,false);
			cache[2*k+1] = calculate_f(w, k,true);
		}
		AliasTable* table = new AliasTable(2*K);
		table->build(cache);
		VecInt& samples = aliasSamples[w];
		VecDouble& props = aliasSampleProps[w];
		if (samples[0] != -1) {
			for (int i = 0; i < nw[w] * ALIAS_SAMPLES_KEPT_RATIO; i++) {
				int sample = table->sample();
				samples[i + 1] = sample;
				props[i + 1] = table->getProportions()[sample];
			}
			samples[0] = 0;
			delete table;
		} else {

			while (aliasTableLocks[w] == 1) {
				printf(
						"WARNING [Table building thread]: Contention while replacing table!\n");
			}
			aliasTableLocks[w] = 1;
			AliasTable* oldTable = aliasTables[w];
			if (oldTable != NULL) {
				delete oldTable;
			}
			aliasTables[w] = table;
			aliasTableLocks[w] = 0;
		}
	}
	virtual void computeTheta(){
		theta.resize(docs.size());
		for(size_t d=0;d<docs.size();d++){
			theta[d].resize(K);
			for(int k=0;k<K;k++){
				theta[d][k]=(getndk(d,k)+alpha)/(docs[d].size()+alphasum);
			}
		}

	}
	virtual void computeThetaSum(){
		for(int k=0;k<K;k++){
			thetasum[k]=0;
		}
		for(size_t d=0;d<docs.size();d++){
			for(int k=0;k<K;k++){
				double val=(getndk(d,k)+alpha)/(docs[d].size()+alphasum);
				thetasum[k]+=val/docs.size();
			}
		}
	}
	virtual void computePhi(){
		phi.resize(V);
		for(int w=0;w<V;w++){
			phi[w].resize(K);
			for(int k=0;k<K;k++){
				double localsum=getphi0(w,k);
				phi[w][k]=(m[w][k]-a*t[w][k])/(b+msum[k])+(a*tsum[k]+b)/(b+msum[k])*localsum;
				phi[w][k]=std::max(0.0, phi[w][k]);
			}
		}
	}
	void update_bucket_dense_sum(int w,int k, bool inc, bool newtable){
		if(inc){
			bucket_dense_sum[w]-=alpha*calculate_f_withoffset(w,k,false,-1,newtable?-1:0);
			bucket_dense_sum[w]-=alpha*calculate_f_withoffset(w,k,true,-1,newtable?-1:0);
			bucket_dense_sum[w]+=alpha*calculate_f(w,k,false);
			bucket_dense_sum[w]+=alpha*calculate_f(w,k,true);
		}else{
			bucket_dense_sum[w]-=alpha*calculate_f_withoffset(w,k,false,1,newtable?1:0);
			bucket_dense_sum[w]-=alpha*calculate_f_withoffset(w,k,true,1,newtable?1:0);
			bucket_dense_sum[w]+=alpha*calculate_f(w,k,false);
			bucket_dense_sum[w]+=alpha*calculate_f(w,k,true);
		}
	}
	int discreteSample(VecInt& indices, VecDouble& prop, int size,
			double roll) {
		for (int i = 0; i < size; i++) {
//			assert(!isnan(prop[i]));
			if (roll <= prop[i])
				return indices[i];
			else
				roll -= prop[i];
		}
		return roll<=MIN_ERROR_BOUND?indices[size-1]:-1;
	}
	void metroHastingSample(int d, int w, int k, int l,bool headOfTable,
			double _bucket_sparse_sum, int discreteSampleResult, int& sample,
			int& sampleType) {
		VecInt& samples = aliasSamples[w];
		VecDouble& sampleProps = aliasSampleProps[w];
		int index = samples[0] + 1;
		int ret = 2*k+headOfTable?1:0;
		double sump = _bucket_sparse_sum + bucket_dense_sum[w];
		bool hasAccept = false;
		while (aliasTableLocks[w] == 1) {
			printf("WARNING: Contention inside core loop\n");
		}
		aliasTableLocks[w] = 1;
		//		while(!aliasTableLocks[w]->try_lock()){printf("WARNING: Contention inside core loop\n");};
		for (int i = 0; i < METRO_HASTING_STEPS; i++) {
			double roll0 = rand01() * sump;
			bool accepted = false;
			if (roll0 <= _bucket_sparse_sum) {
				sample = discreteSampleResult;
				assert(sample != -1);
				sampleType = 0;
				accepted = true;
//								printf("*MH=%d d=%d l=%d w=%d k=%d sample=%d "
//										"type=%d sparse_sum=%lf dense_sum=%lf roll=%lf\n",i,d,l,w,k,sample,
//										sampleType,_bucket_sparse_sum,bucket_dense_sum[w],roll0);
			} else {
				double newZSampleProp;
				if (index != 0) {
					sample = samples[index];
					newZSampleProp = sampleProps[index];
					index++;
				} else {
					AliasTable* table = aliasTables[w];
					assert(table!=NULL);
					sample = table->sample();
					newZSampleProp = table->getProportions()[sample];
				}
				sampleType = 1;
				double propPrevSample = calculate_f_flat(w, ret);
				double propNewSample = calculate_f_flat(w, sample);
				int ndk0 = getndk(d, ret/2);
				int ndk1 = getndk(d, sample/2);
				double acceptance =
						((ndk1 + alpha) / (ndk0 + alpha))
								* (propNewSample / propPrevSample)
								* ((alpha * propPrevSample
										+ propPrevSample * ndk0)
										/ (alpha * newZSampleProp
												+ propNewSample * ndk1));
				double roll = rand01();
				if (roll < acceptance) {
					accepted = true;
				}else{
					numMetroRej++;
				}
//								printf("MH=%d d=%d l=%d w=%d oldsample=%d sample=%d "
//										"type=%d sparse_sum=%lf dense_sum=%lf ndk0=%d "
//										"ndk1=%d newZSampleProp=%lf propPrevSample=%lf propNewSample=%lf acceptance=%.13lf roll=%lf accepted=%d\n",
//										i,d,l,w,ret,sample,
//										sampleType, _bucket_sparse_sum,bucket_dense_sum[w], ndk0,
//										ndk1,newZSampleProp,propPrevSample,propNewSample,acceptance,roll,accepted);

			}
			if (accepted) {
				ret = sample;
				hasAccept = true;
			}
			assert(index <= nw[w] * ALIAS_SAMPLES_KEPT_RATIO);
		}
		if (index != 0)
			samples[0] = index - 1;
		aliasTableLocks[w] = 0;
		//		aliasTableLocks[w]->unlock();
		if (!hasAccept)
			sampleType = 2;

		sample = ret;
	}
	virtual void sampleWord(int d, int l) {
		int w = docs[d][l];
		int k = z[d][l];
		assert(m[w][k]>=t[w][k]);
		assert(!(t[w][k]==0 && m[w][k]!=0));
		bool removeTable=(rand01()<=((double)(t[w][k])/(double)m[w][k])?true:false);
		dec_m(w,k);
		dec_ndk(d,k);
		removeTable=removeTable&& (!(t[w][k]==1 && m[w][k]!=0));
		if(removeTable){
			dec_t(w,k);
		}
		if(t[w][k]==0 && m[w][k]!=0){
			printf("error: d=%d l=%d t[%d][%d]=%d m[%d][%d]=%d\n",d,l,w,k,t[w][k],w,k,m[w][k]);
		}
		assert(!(t[w][k]==0 && m[w][k]!=0));
		update_bucket_dense_sum(w, k, false,removeTable);
		double bucket_sparse_sum = 0;
		VecInt& ndkSparseLocal = ndkSparse[d];
		for (int i = 0; i < ndkSparseLocal.size(); i++) {
			int encodedValue = ndkSparseLocal[i];
			int topic = decodeTopic(encodedValue);
			int count = decodeCount(encodedValue);
			double val0 = count * calculate_f(w,topic,false);
			bucket_sparse_values[2*i] = val0;
			bucket_sparse_indices[2*i] = 2*topic;
			bucket_sparse_sum += val0;
			double val1 = count * calculate_f(w,topic,true);
			bucket_sparse_values[2*i+1] = val1;
			bucket_sparse_indices[2*i+1] = 2*topic+1;
			bucket_sparse_sum += val1;

		}
		double roll=rand01() * bucket_sparse_sum;
		int discreteSampleResult = discreteSample(bucket_sparse_indices,
				bucket_sparse_values, ndkSparseLocal.size()*2,
				roll);
		//		sparsePartTable->build(bucket_sparse_values);
		assert(discreteSampleResult>=0);
		assert(discreteSampleResult<2*K);
		if(discreteSampleResult<=-1){
			printf(" d=%d,l=%d | oldtopic=%d newsample=%d sparse_sum=%lf roll=%lf\n",d,l ,k, discreteSampleResult,  bucket_sparse_sum,roll );
			for(int i=0;i<ndkSparseLocal.size()*2;i++){
				printf("--bucket[%d]=%lf topic=%d r=%d, fval=%lf\n",i,bucket_sparse_values[i], bucket_sparse_indices[i]/2,bucket_sparse_indices[i]%2,calculate_f(w,i/2,i%2==1));
			}
		}
		int newSample = -1, flag = -1;
		metroHastingSample(d, w, k, l,removeTable,
				bucket_sparse_sum, discreteSampleResult, newSample,
				flag);
		sampleFlagCount[flag]++;
		//		if(newTopic==-1 || flag==-1){
		//			printf(" d=%d,l=%d | flag=%d oldtopic=%d newtopic=%d\n",d,l,flag ,k, newTopic);
		//		}
		assert(newSample >= 0);
		assert(newSample < 2*K);

		int newTopic=newSample/2;
		int headOfTable=newSample%2;
		assert(newTopic<K);
		inc_m(w,newTopic);
		inc_ndk(d, newTopic);
		if(headOfTable==1){
			inc_t(w, newTopic);
		}
		update_bucket_dense_sum(w, newTopic, true,headOfTable);

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
			printf("flags %d=%d\n",i,sampleFlagCount[i]);
			sampleFlagCount[i]=0;
		}
		printf("average numMetroRej per Metro step=%lf\n",numMetroRej/(double)METRO_HASTING_STEPS);
		numMetroRej=0;
	}

	virtual void setAlphaVec(const VecDouble& alphaVec) {
		this->alphaVec = alphaVec;
	}
	virtual void setAlpha(double alpha) {
		this->alpha = alpha;
	}
	virtual double getGamma() const {
		return gamma;
	}
	virtual void setGamma(double gamma) {
		this->gamma = gamma;
	}
	virtual void setDiscount(double a) {
		this->a = a;
	}
	virtual double getDiscount() {
		return a;
	}
	virtual void setConcentration(double b) {
		this->b = b;
	}
	virtual double getConcentration() {
		return b;
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

	virtual void clearPhi() {
		phi.clear();
	}
	virtual void clearTheta() {
		theta.clear();
	}
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
	virtual const VecDouble& getAlphaVec() const {
		return alphaVec;
	}
	virtual void copyState(ILDA* _lda){
		AliasWordPDPLDA* lda=static_cast<AliasWordPDPLDA*> (_lda);
		K=lda->K;
		V=lda->V;
		m=lda->m;
		msum=lda->msum;
		t=lda->t;
		tsum=lda->tsum;
		a=lda->a;
		b=lda->b;
		gamma=lda->gamma;
		gammasum=lda->gammasum;
		gsnCache=lda->gsnCache;
		gsn_maxM=lda->gsn_maxM;
		gsn_maxN=lda->gsn_maxN;
		alpha=lda->alpha;
		alphasum=lda->alphasum;
		alphaVec=lda->alphaVec;
	}

};
}

