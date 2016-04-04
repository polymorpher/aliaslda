#pragma once
#include "GeneralizedStirlingNumber.hpp"
#include "ILDA.hpp"
namespace AliasLDA {
class AliasHDPLDA: public ILDA {
	friend class HDPLDA;
	const double MIN_ERROR=1e-9;
	Vector2DInt docs;
	int V;
	int K;
	double b0; //concentration
	double b1;
	double gamma;  //for generating base Dirichlet
	double gammasum;
	Vector2DDouble phi; //dim K V. Word|topic,group distribution.
	Vector2DDouble theta; //dim I D K topic|document,group distribution.

	VecDouble thetasum; //dim I K theta_{I,D,K} sum over all documents

	Vector2DInt mwkSparse; //dim V K. ==nwk;
	Vector2DInt mwkIndex;
	VecInt mk; //dim K.
	int m0;
	Vector2DInt tdk;  //dim V K. multiplicity
	VecInt tk; //dim K
	int tsum;
	VecInt t0k;
	Vector2DInt ndk;
	VecInt n0k;
	Vector2DInt z; //dim I _D __L. topic assignment
	Vector2DInt u;
	VecInt nw;
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;
	GeneralizedStirlingNumber& gsn;
	Vector2DDouble gsnCache;
	int gsn_maxN;
	int gsn_maxM;

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
	Vector2DDouble bucket_dense_values;
	//threads
	std::thread* thread_table;
	bool terminated;
	VecInt aliasTableLocks;
	std::vector<AliasTable*> aliasTables;
	const int METRO_HASTING_STEPS = 2;
	const int ALIAS_SAMPLES_KEPT_RATIO = 16;
	const bool USE_MIXED_SUPPLY_CONSUME_MECH=true;
	int numMetroRej;
	void initialiseGSN(){
		NMPair limits=gsn.getLimit(0);
		while(limits.m==-1)limits=gsn.getLimit(0);
		gsn_maxN=limits.n;
		gsn_maxM=limits.m;
		if(limits.m<gsn_maxM||limits.n<gsn_maxN){
			printf("initialising GSN Module...\n");
			while(!gsn.initialize(0,gsn_maxN,gsn_maxM)){}
			printf("GSN Module initialized\n");
		}
		double* tempCache;
		gsnCache.resize(gsn_maxN);
		while((tempCache=gsn.getCache(0))==NULL){}
		printf("Copying over GSN to local buffer\n");
		for(int i=0;i<limits.n;i++){
			gsnCache[i].resize(gsn_maxM);
			for(int j=0;j<std::min(limits.m,i);j++){
				int ind=i * limits.m + j;
				gsnCache[i][j]=tempCache[ind];
//				printf("gsn %d %d val=%lf\n",i,j,gsnCache[i][j]);
			}
		}
	}
	double getGSNRatio(int n0,int m0,int n1,int m1){
		assert(n0<gsn_maxN && n1<gsn_maxN);
		assert(m0<gsn_maxM && m1<gsn_maxM);
		assert(gsnCache[n0][m0]!=NAN && gsnCache[n1][m1]!=NAN);
//		if(gsnCache[n0][m0] < gsnCache[n1][m1] ){
//		printf("lS(%d,%d)=%lf, lS(%d,%d)=%lf, ratio=%lf\n",n0,m0,gsnCache[n0][m0],n1,m1,gsnCache[n1][m1],exp(gsnCache[n0][m0]-gsnCache[n1][m1]));
//		}
		return exp(gsnCache[n0][m0]-gsnCache[n1][m1]);
	}
public:
	inline int getmwk(int w,int k){
		int index=mwkIndex[w][k];
		return index==0?0:decodeCount(mwkSparse[w][index-1]);
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

	inline void inc_mwk(int w, int k) {
		int index = mwkIndex[w][k] - 1;
		if (index < 0) {
			mwkSparse[w].push_back(oneWithTopicBitShift + k);
			index = mwkSparse[w].size() - 1;
			mwkIndex[w][k] = index + 1;
		} else {
			incEncodedTopic (mwkSparse[w][index]);
		}
		oneStepSortInc(mwkIndex[w], mwkSparse[w], index);
		mk[k]++;
	}
	inline bool dec_mwk(int w, int k) {
		int index = mwkIndex[w][k] - 1;
		assert(index >= 0);
		assert(mwkSparse[w][index] >= 0);
		decEncodedTopic (mwkSparse[w][index]);
		assert(mwkSparse[w].size() > index);
		mk[k]--;
		return oneStepSortDec(mwkIndex[w], mwkSparse[w], index);
	}

	inline void inc_tdk(int d, int k) {
		tdk[d][k]++;
		tk[k]++;
		tsum++;
	}
	inline void dec_tdk(int d, int k) {
		assert(tdk[d][k]>0);
		tdk[d][k]--;
		tk[k]--;
		tsum--;
	}
	inline void inc_t0k(int k) {
		t0k[k]++;
		tk[k]++;
		tsum++;
	}
	inline void dec_t0k( int k) {
		t0k[k]--;
		tk[k]--;
		tsum--;
	}
	inline void inc_n0k(int k){
		n0k[k]++;
	}
	inline void dec_n0k(int k){
		n0k[k]--;
	}
	inline void inc_ndk(int d, int k) {
		ndk[d][k]++;
	}
	inline void dec_ndk(int d, int k) {
		assert(ndk[d][k]>0);
		ndk[d][k]--;
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
	AliasHDPLDA() :
			K(0), V(0), u01(0, 1), gsn(GeneralizedStirlingNumber::getInstance()), terminated(false) {
		rgen.seed(time(NULL) + rand());
	}
	virtual ~AliasHDPLDA(){
		terminated=true;
		thread_table->join();
	}
	inline double rand01() {
		return u01(rgen);
	}
	inline int randInt(int limit) {
		return (int) (rand01() * limit);
	}


	void initialiseWordCount() {
		for (int d = 0; d < docs.size(); d++) {
			for (int l = 0; l < docs[d].size(); l++) {
				nw[docs[d][l]]++;
			}
		}
	}
	void initialiseWordRestaurantCounter() {
		mwkSparse.resize(V);
		mwkIndex.resize(V);
		mk.resize(K);
		m0=0;
		t0k.resize(K);
		n0k.resize(K);
		for (int w = 0; w < V; w++) {
			mwkSparse[w].reserve(std::min(K,nw[w]));
			mwkIndex[w].resize(K);
		}
	}
	void initialiseTopicRestaurantCounter(){
		tdk.resize(docs.size());
		tk.resize(K);
		for (int d = 0; d < docs.size(); d++) {

			tdk[d].resize(K);
		}
	}

	void initialiseDocumentTopicCounter() {
		ndk.resize(docs.size());
		for (int d = 0; d < docs.size(); d++) {
			ndk[d].resize(K);
		}
	}
	 double calculate_f_all_in_one_flat(int d,int composite){
		return calculate_f_all_in_one(d,composite/2,composite%2==1);
	}
	 double calculate_f_all_in_one(int d, int k, bool u){
		if(u){
			if(t0k[k]==0){
//				assert(mk[k]==0);
				return (b0*b1)/(b0+tsum)/(gammasum+mk[k]);

			}else{
				if(tdk[d][k]==0){
					return (b1*tk[k])/(double)(tk[k]+1)*((tk[k])/(tsum+b0))/(gammasum+mk[k]);
				}else{
					return getGSNRatio(ndk[d][k]+1,tdk[d][k]+1,ndk[d][k],tdk[d][k])*
							(tdk[d][k]+1)/(ndk[d][k]+1)/
							(gammasum+mk[k]);
				}
			}

		}else{
			if(t0k[k]!=0){
				return 0;
			}else{
				return (b0*b1)/(b0+tsum)/(gammasum+mk[k]);
			}
		}
	}
	 double calculate_f_all_in_one_with_offset(int d, int k, bool u, int noffset,int toffset,int t0offset, int moffset){
		if(u){
			if(t0k[k]+t0offset==0){
				return (b0*b1)/(b0+tsum+toffset+t0offset)/(gammasum+mk[k]+moffset);
			}else{
				if(tdk[d][k]+toffset==0){
					return (b1*(tk[k]+toffset+t0offset))/(double)(tk[k]+toffset+t0offset+1)*((tk[k]+toffset+t0offset)/(tsum+toffset+t0offset+b0))/(gammasum+mk[k]+moffset);
				}else{
					return getGSNRatio(ndk[d][k]+1+noffset,tdk[d][k]+1+toffset,ndk[d][k]+noffset,tdk[d][k]+toffset)*
							(tdk[d][k]+1+toffset)/(ndk[d][k]+noffset+1)/
							(gammasum+mk[k]+moffset);
				}
			}
		}else{
			if(t0k[k]+t0offset!=0){
				return 0;
			}else{
				return (b0*b1)/(b0+tsum+toffset+t0offset)/(gammasum+mk[k]+moffset);
			}
		}
	}
	void initialiseGamma(){
		gammasum=gamma*V;
	}
	void initialiseTopicsAndCounts() {
		z.resize(docs.size());
		u.resize(docs.size());
		for (size_t d = 0; d < docs.size(); d++) {
			VecInt& doc = docs[d];
			z[d].resize(doc.size());
			u[d].resize(doc.size());
			VecInt& zd = z[d];
			VecInt& ud = u[d];
			for (size_t l = 0; l < doc.size(); l++) {
				int w = doc[l];
				int k = randInt(K);
				inc_ndk(d, k);
				inc_mwk(w, k);
				zd[l] = k;
				ud[l]=1;
				if(tdk[d][k]==0){
					inc_tdk(d,k);
					ud[l]=1;
				}else{
					if(rand01() <= 1/(double)(tdk[d][k]+1)){
						inc_tdk(d,k);
						if(rand01() <= 1/(double)(t0k[k]+1)){
							inc_t0k(k);
						}
						inc_n0k(k);
						m0++;
						ud[l]=1;
					}
				}
				if(t0k[k]==0){
					inc_t0k(k);
					inc_n0k(k);
					m0++;
					ud[l]=0;
				}


				assert(!(t0k[k]==0 && n0k[k]!=0));
				assert(!(tdk[d][k]==0 && getmwk(w,k)!=0));
			}
		}
		initialise_bucket_dense_sum();
	}
	void initialise_bucket_dense_sum(){
		bucket_dense_values.resize(docs.size());
		bucket_dense_sum.resize(docs.size());
		for(int d=0;d<docs.size();d++){
			bucket_dense_sum[d]=0;
			bucket_dense_values[d].resize(2*K);
			for(int k=0;k<K;k++){
				double val0=gamma*calculate_f_all_in_one(d,k,false);
				double val1=gamma*calculate_f_all_in_one(d,k,true);
				bucket_dense_values[d][2*k]=val0;
				bucket_dense_values[d][2*k+1]=val1;
				bucket_dense_sum[d]+=val0;
				bucket_dense_sum[d]+=val1;

			}
		}
	}
	void initialise_bucket_sparse(){
		bucket_sparse_values.resize(2*K);
		bucket_sparse_indices.resize(2*K);
	}
	void initialise_sampleFlagCount(){
		sampleFlagCount[0]=0;
		sampleFlagCount[1]=0;
		sampleFlagCount[2]=0;
	}
	void initialiseAliasTables() {
		aliasTables.resize(docs.size());
		aliasSamples.resize(docs.size());
		aliasSampleProps.resize(docs.size());
		aliasTableLocks.resize(docs.size());
		for (int d = 0; d < aliasSamples.size(); d++) {
			if (USE_MIXED_SUPPLY_CONSUME_MECH
					&& (docs[d].size() * ALIAS_SAMPLES_KEPT_RATIO
							> sizeof(AliasTable)
									+ sizeof(AliasTable::AliasElement) * K)) {
				aliasSamples[d].resize(1);
				aliasSamples.shrink_to_fit();
				aliasSampleProps[d].resize(1);
				aliasSampleProps.shrink_to_fit();
				aliasSamples[d][0] = -1;
				aliasSampleProps[d][0] = -1;
			} else {
				aliasSamples[d].resize(docs[d].size() * ALIAS_SAMPLES_KEPT_RATIO + 1);
				aliasSampleProps[d].resize(
						docs[d].size() * ALIAS_SAMPLES_KEPT_RATIO + 1);
			}
		}
		for (int d = 0; d < docs.size(); d++)
			computeAliasTable(d);
		thread_table = new std::thread(
				&AliasLDA::AliasHDPLDA::thread_computeAliasTable, this);
	}
	virtual void initialise() {
		initialiseGSN();
		initialiseGamma();
		nw.resize(V);
		initialiseWordCount();
		initialiseWordRestaurantCounter();
		initialiseAssignmentOnly();
		printf("initialised.\m");
	}
	virtual void initialiseAssignmentOnly(){
		thetasum.resize(K);
		initialiseDocumentTopicCounter();
		initialiseTopicRestaurantCounter();
		initialiseBitShifts();
		initialise_bucket_sparse();
		initialise_sampleFlagCount();
		initialiseTopicsAndCounts();
		initialiseAliasTables();
	}
	void computeAliasTable(int d){
		VecDouble cache;
		cache.resize(2*K);
		for (int k = 0; k < K; k++) {
			cache[2*k] = calculate_f_all_in_one(d, k,false);
			cache[2*k+1] = calculate_f_all_in_one(d, k,true);
		}
		AliasTable* table = new AliasTable(2*K);
		table->build(cache);
		//if(nw[w]*ALIAS_SAMPLES_KEPT_RATIO < sizeof(AliasTable)+sizeof(AliasTable::AliasElement)*K){

		VecInt& samples=aliasSamples[d];
		VecDouble& props=aliasSampleProps[d];
		if(samples[0]!=-1){
			for(int i=0;i<docs[d].size()*ALIAS_SAMPLES_KEPT_RATIO;i++){
				int sample=table->sample();
				samples[i+1]=sample;
				props[i+1]=table->getProportions()[sample];
			}
			samples[0]=0;
			delete table;
		}else{
//			std::shared_ptr<AliasTable> oldTable=aliasTables[w];
			while(aliasTableLocks[d]==1){printf("WARNING [Table building thread]: Contention while replacing table!\n");}
			aliasTableLocks[d]=1;
//			while(!aliasTableLocks[w]->try_lock()){printf("WARNING [Table building thread]: Contention while replacing table!\n");}
			AliasTable* oldTable=aliasTables[d];
			if(oldTable!=NULL){
				delete oldTable;
			}
			aliasTables[d]=table;
//			aliasTableLocks[w]->unlock();
			aliasTableLocks[d]=0;
		}
	}
	void thread_computeAliasTable(){
		while(!terminated){
//			printf("rebuilding aliastable...\n");
			for(int d=0;d<docs.size();d++){
				computeAliasTable(d);
			}
//			printf("rebuilding aliastable complete.\n");
		}
	}
	inline double gettheta0(int d,int k){
		return n0k[k]/(b0+m0)+b0/(b0+m0)*(1/K);
	}
	virtual void computeTheta(){
		theta.resize(docs.size());
		for(size_t d=0;d<docs.size();d++){
			theta[d].resize(K);
			for(int k=0;k<K;k++){
				theta[d][k]=std::max(0.0,ndk[d][k]*tdk[d][k]/(b1+docs[d].size())+(b1)/(b1+docs[d].size())*gettheta0(d,k));
			}
		}

	}
	virtual void computeThetaSum(){
		for(int k=0;k<K;k++){
			thetasum[k]=0;
		}
		for(size_t d=0;d<docs.size();d++){
			for(int k=0;k<K;k++){
				double val=ndk[d][k]/(b1+docs[d].size())+(b1)/(b1+docs[d].size())*gettheta0(d,k);
				thetasum[k]+=std::max(0.0,val)/docs.size();
			}
		}
	}

	virtual void computePhi(){
		phi.resize(V);
		for(int w=0;w<V;w++){
			phi[w].resize(K);
			for(int k=0;k<K;k++){
				phi[w][k]=(gamma+getmwk(w,k))/(gammasum+mk[k]);
				phi[w][k]=std::max(0.0, phi[w][k]);
			}
		}
	}

	int discreteSample(VecInt& indices, VecDouble& prop, int size,
			double roll) {
		for (int i = 0; i < size; i++) {
			if (roll <= prop[i])
				return indices[i];
			else
				roll -= prop[i];
		}
		if(roll<MIN_ERROR)return size-1;
		return -1;
	}
	void metroHastingSample(int d, int w, int k, int l, int ul,
			double _bucket_sparse_sum, int discreteSampleResult, int& sample,
			int& sampleType) {
		VecInt& samples = aliasSamples[d];
		VecDouble& sampleProps = aliasSampleProps[d];
		int index = samples[0] + 1;
		int ret = 2 * k + ul;
		double sump = _bucket_sparse_sum + bucket_dense_sum[d];
		bool hasAccept = false;
		while (aliasTableLocks[d] == 1) {
			printf("WARNING: Contention inside core loop\n");
		}
		aliasTableLocks[d] = 1;
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
					AliasTable* table = aliasTables[d];
					assert(table!=NULL);
					sample = table->sample();
					newZSampleProp = table->getProportions()[sample];
				}
				sampleType = 1;
				double propPrevSample = calculate_f_all_in_one_flat(d, ret);
				double propNewSample = calculate_f_all_in_one_flat(d, sample);
				int mwk0 = getmwk(w, ret / 2);
				int mwk1 = getmwk(w, sample / 2);
				double acceptance =
						((mwk1 + gamma) / (mwk0 + gamma))
								* (propNewSample / propPrevSample)
								* ((gamma * propPrevSample
										+ propPrevSample * mwk0)
										/ (gamma * newZSampleProp
												+ propNewSample * mwk1));
				double roll = rand01();
				if (roll < acceptance) {
					accepted = true;
				} else {
					numMetroRej++;
				}
//				printf("MH=%d d=%d l=%d w=%d oldsample=%d sample=%d "
//						"type=%d sparse_sum=%lf dense_sum=%lf ndk0=%d "
//						"ndk1=%d newZSampleProp=%lf propPrevSample=%lf propNewSample=%lf acceptance=%.13lf roll=%lf accepted=%d\n",
//						i,d,l,w,ret,sample,
//						sampleType, _bucket_sparse_sum,bucket_dense_sum[w], mwk0,
//						mwk1,newZSampleProp,propPrevSample,propNewSample,acceptance,roll,accepted);

			}
			if (accepted) {
				ret = sample;
				hasAccept = true;
			}
			assert(index <= docs[d].size() * ALIAS_SAMPLES_KEPT_RATIO);
		}
		if (index != 0)
			samples[0] = index - 1;
		aliasTableLocks[d] = 0;
		//		aliasTableLocks[w]->unlock();
		if (!hasAccept)
			sampleType = 2;

		sample = ret;
	}
	virtual void sampleWord(int d, int l) {
		int w = docs[d][l];
		int k = z[d][l];
		int ul = u[d][l];
		bool removeRootTable=false;
		bool removeTable=false;
		if(ul==0){
//			removeRootTable = (rand01()<=((double)(t0k[k])/(double)n0k[k])?true:false);
			removeRootTable=true;
			dec_n0k(k);
			m0--;
		}
//		removeTable = (rand01()<=((double)(tdk[d][k])/(double)(ndk[d][k]))?true:false);
		removeTable=true;
		dec_ndk(d,k);
		dec_mwk(w,k);
		removeRootTable=removeRootTable&&(!(t0k[k]==1 && n0k[k]!=0)) && t0k[k]>0;
		removeTable=removeTable&&(!(tdk[d][k]==1 && ndk[d][k]!=0))  && tdk[d][k]>0;
//		printf("error: d=%d l=%d k=%d ndk[%d][%d]=%d tdk[%d][%d]=%d rmTable=%d rmRoot=%d\n",d,l,k,d,k,ndk[d][k],d,k,tdk[d][k], removeTable, removeRootTable);

		if(removeRootTable)dec_t0k(k);
		if(removeTable)dec_tdk(d,k);
//		if(ndk[d][k]<tdk[d][k]){
//			printf("error: d=%d l=%d k=%d ndk[%d][%d]=%d tdk[%d][%d]=%d rmTable=%d rmRoot=%d\n",d,l,k,d,k,ndk[d][k],d,k,tdk[d][k], removeTable, removeRootTable);
//		}
		assert(n0k[k]>=t0k[k]);
		assert(!(t0k[k]==0 && n0k[k]!=0));
		assert(ndk[d][k]>=tdk[d][k]);
		if((tdk[d][k]==0 && ndk[d][k]!=0))
			printf("error: d=%d l=%d k=%d ndk[%d][%d]=%d tdk[%d][%d]=%d rmTable=%d rmRoot=%d\n",d,l,k,d,k,ndk[d][k],d,k,tdk[d][k], removeTable, removeRootTable);
		assert(!(tdk[d][k]==0 && ndk[d][k]!=0));
//
//		if(t[w][k]==0 && m[w][k]!=0){
//			printf("error: d=%d l=%d t[%d][%d]=%d m[%d][%d]=%d\n",d,l,w,k,t[w][k],w,k,m[w][k]);
//		}
		update_bucket_dense_sum(d, k, false ,removeRootTable,removeTable);
		double bucket_sparse_sum = 0;
		VecInt& mwkSparseLocal = mwkSparse[w];
		for (int i = 0; i < mwkSparseLocal.size(); i++) {
			int encodedValue = mwkSparseLocal[i];
			int topic = decodeTopic(encodedValue);
			int count = decodeCount(encodedValue);
			double val0 = count * calculate_f_all_in_one(d,topic,false);
			bucket_sparse_values[2*i] = val0;
			bucket_sparse_indices[2*i] = 2*topic;
			bucket_sparse_sum += val0;
			double val1 = count * calculate_f_all_in_one(d,topic,true);
			bucket_sparse_values[2*i+1] = val1;
			bucket_sparse_indices[2*i+1] = 2*topic+1;
			bucket_sparse_sum += val1;

		}
		int discreteSampleResult=-1;
		if(mwkSparseLocal.size()>0){
			double roll=rand01() * bucket_sparse_sum;
			discreteSampleResult= discreteSample(bucket_sparse_indices,
					bucket_sparse_values, mwkSparseLocal.size()*2,
					roll);
			//		sparsePartTable->build(bucket_sparse_values);

			if(discreteSampleResult < 0){
				printf(" d=%d,l=%d | oldtopic=%d newsample=%d sparse_sum=%lf roll=%lf\n",d,l ,k, discreteSampleResult,  bucket_sparse_sum,roll );
				for(int i=0;i<mwkSparseLocal.size()*2;i++){
					printf("--bucket[%d]=%lf topic=%d r=%d, fval=%lf\n",i,bucket_sparse_values[i], bucket_sparse_indices[i]/2,bucket_sparse_indices[i]%2,calculate_f_all_in_one(d,i/2,i%2==1));
				}
			}
			assert(discreteSampleResult>=0);
			assert(discreteSampleResult<2*K);
		}
		int newSample = -1, flag = -1;
		metroHastingSample(d, w, k, l,ul,
				bucket_sparse_sum, discreteSampleResult, newSample,
				flag);
		sampleFlagCount[flag]++;
		//		if(newTopic==-1 || flag==-1){
		//			printf(" d=%d,l=%d | flag=%d oldtopic=%d newtopic=%d\n",d,l,flag ,k, newTopic);
		//		}
		assert(newSample >= 0);
		assert(newSample < 2*K);

		int newTopic=newSample/2;
		int newul=newSample%2;
		inc_mwk(w,newTopic);
		inc_ndk(d, newTopic);

		bool addRootTable=false,addTable=false;
//		addRootTable=(rand01()<=(1-(double)(t0k[newTopic])/(double)n0k[newTopic])?true:false);
		addRootTable=(addRootTable || t0k[newTopic]==0) && newul==0;
//		addTable=(rand01()<=(1-(double)(tdk[d][newTopic])/(double)ndk[d][newTopic])?true:false);
		addTable= (addTable || tdk[d][newTopic]==0) ;
		if(newul==0){
			inc_n0k(newTopic);
			m0++;
		}
		if(addRootTable){inc_t0k(newTopic);}
		if(addTable){
			inc_tdk(d, newTopic);
		}
		update_bucket_dense_sum(d, newTopic, true ,addRootTable,addTable);
		u[d][l] = newul;
		z[d][l] = newTopic;
	}
//	void update_bucket_dense_sum(int d,int k,  bool inc, bool newRootTable,bool newTable){
//		if(inc){
//			bucket_dense_sum[d]-=gamma*calculate_f_all_in_one_with_offset(d,k,false,-1,newTable?-1:0,newRootTable?-1:0,-1);
//			bucket_dense_sum[d]-=gamma*calculate_f_all_in_one_with_offset(d,k,true,-1,newTable?-1:0,newRootTable?-1:0,-1);
//			bucket_dense_sum[d]+=gamma*calculate_f_all_in_one(d,k,false);
//			bucket_dense_sum[d]+=gamma*calculate_f_all_in_one(d,k,true);
//		}else{
//			bucket_dense_sum[d]-=gamma*calculate_f_all_in_one_with_offset(d,k,false,1, newTable?1:0, newRootTable?1:0 ,1);
//			bucket_dense_sum[d]-=gamma*calculate_f_all_in_one_with_offset(d,k,true,1, newTable?1:0, newRootTable?1:0 ,1);
//			bucket_dense_sum[d]+=gamma*calculate_f_all_in_one(d,k,false);
//			bucket_dense_sum[d]+=gamma*calculate_f_all_in_one(d,k,true);
//		}
//	}
	void update_bucket_dense_sum(int d,int k,  bool inc, bool newRootTable,bool newTable){
		bucket_dense_sum[d]-=bucket_dense_values[d][2*k];
		bucket_dense_sum[d]-=bucket_dense_values[d][2*k+1];
		double newval0=calculate_f_all_in_one(d,k,false);
		double newval1=gamma*calculate_f_all_in_one(d,k,true);
		bucket_dense_sum[d]+=newval0;
		bucket_dense_sum[d]+=newval1;
		bucket_dense_values[d][2*k]=newval0;
		bucket_dense_values[d][2*k+1]=newval1;
	}
	virtual void gibbsStep() {
		numMetroRej=0;
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
	}


	virtual double getGamma() const {
		return gamma;
	}
	virtual void setGamma(double gamma) {
		this->gamma = gamma;
	}
	virtual void setConcentration(double b) {
		this->b1 = b;
	}
	virtual double getConcentration() {
		return b1;
	}
	virtual void setRootConcentration(double b) {
		this->b0 = b;
	}
	virtual double getRootConcentration() {
		return b0;
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
	virtual void copyState(ILDA* _lda){
		AliasHDPLDA* lda=static_cast<AliasHDPLDA*> (_lda);
		K=lda->K;
		V=lda->V;
		mwkSparse=lda->mwkSparse;
		mwkIndex=lda->mwkIndex;
		mk=lda->mk;
		n0k=lda->n0k;
		t0k=lda->t0k;
		tk=lda->tk;
		tsum=lda->tsum;
		m0=lda->m0;
		b0=lda->b0;
		b1=lda->b1;
		gamma=lda->gamma;
		gammasum=lda->gammasum;
		gsnCache=lda->gsnCache;
		gsn_maxM=lda->gsn_maxM;
		gsn_maxN=lda->gsn_maxN;
	}
};
}

