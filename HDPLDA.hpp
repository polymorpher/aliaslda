#pragma once
#include "GeneralizedStirlingNumber.hpp"
#include "ILDA.hpp"
#include"AliasHDPLDA.hpp"
namespace AliasLDA {
class HDPLDA: public ILDA {

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

	Vector2DInt mwk; //dim V K. ==nwk;
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
	VecDouble bucket;
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;
	GeneralizedStirlingNumber& gsn;
	Vector2DDouble gsnCache;
	int gsn_maxN;
	int gsn_maxM;

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


	inline void inc_mwk(int w, int k) {	//increment m
		mwk[w][k]++;
		mk[k]++;
	}
	inline void dec_mwk(int w, int k) {	//decrement m
		mwk[w][k]--;
		mk[k]--;
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

public:
	HDPLDA() :
			K(0), V(0), u01(0, 1), gsn(GeneralizedStirlingNumber::getInstance()) {
		rgen.seed(time(NULL) + rand());
	}
	inline double rand01() {
		return u01(rgen);
	}
	inline int randInt(int limit) {
		return (int) (rand01() * limit);
	}

	inline int getndk(int d, int k) {
		return ndk[d][k];
	}
	void initialiseWordCount() {
		for (int d = 0; d < docs.size(); d++) {
			for (int l = 0; l < docs[d].size(); l++) {
				nw[docs[d][l]]++;
			}
		}
	}
	void initialiseWordRestaurantCounter() {
		mwk.resize(V);
		mk.resize(K);
		m0=0;
		t0k.resize(K);
		n0k.resize(K);
		for (int w = 0; w < V; w++) {
			mwk[w].resize(K);
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
	inline double calculate_f_all_in_one_flat(int w,int composite){
		return calculate_f_all_in_one(w,composite/2,composite%2==1);
	}
	inline double calculate_f_all_in_one(int d, int k, bool u){
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
	inline double calculate_f_all_in_one_with_offset(int d, int k, bool u, int noffset,int toffset, int moffset){
		if(u){
			if(t0k[k]==0){
				return 0;
			}else{
				if(tdk[d][k]==0){
					return (b1*(tk[k]+toffset))/(double)(tk[k]+toffset+1)*((tk[k]+toffset)/(tsum+toffset+b0))/(gammasum+mk[k]+moffset);
				}else{
					return getGSNRatio(ndk[d][k]+1+noffset,tdk[d][k]+1+toffset,ndk[d][k]+noffset,tdk[d][k]+toffset)*
							(tdk[d][k]+1+toffset)/(ndk[d][k]+noffset+1)/
							(gammasum+mk[k]+moffset);
				}
			}

		}else{
			if(t0k[k]!=0){
				return 0;
			}else{
				return (b0*b1)/(b0+tsum)/(gammasum+mk[k]+moffset);
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
				assert(!(tdk[d][k]==0 && mwk[w][k]!=0));
			}
		}
	}
	virtual void initialise() {
		initialiseGSN();
		initialiseGamma();
		initialiseWordRestaurantCounter();
		initialiseAssignmentOnly();
		printf("initialised.\m");
	}
	virtual void initialiseAssignmentOnly(){
		thetasum.resize(K);
		nw.resize(V);
		initialiseDocumentTopicCounter();
		initialiseWordCount();
		initialiseTopicRestaurantCounter();
		bucket.resize(2*K);
		initialiseTopicsAndCounts();
	}
	inline double gettheta0(int d,int k){
		return n0k[k]/(b0+m0)+b0/(b0+m0)*(1/K);
	}
	virtual void computeTheta(){
		theta.resize(docs.size());
		for(size_t d=0;d<docs.size();d++){
			theta[d].resize(K);
			for(int k=0;k<K;k++){
				theta[d][k]=std::max(0.0,ndk[d][k]/(b1+docs[d].size())+(b1)/(b1+docs[d].size())*gettheta0(d,k));
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
				phi[w][k]=(gamma+mwk[w][k])/(gammasum+mk[k]);
				phi[w][k]=std::max(0.0, phi[w][k]);
			}
		}
	}


	int discreteSample( VecDouble& prop, int size,
			double roll) {
		for (int i = 0; i < size; i++) {
			if (roll <= prop[i])
				return i;
			else
				roll -= prop[i];
		}
		return -1;
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
		double bucket_sparse_sum = 0;
		VecInt& mwkLocal = mwk[w];
		for (int i = 0; i < K; i++) {
			double val0 = (mwkLocal[i] +gamma) * calculate_f_all_in_one(d,i,false);
			bucket[2*i] = val0;
			bucket_sparse_sum += val0;

			double val1 = (mwkLocal[i]+ gamma)  * calculate_f_all_in_one(d,i,true);
			bucket[2*i+1] = val1;
			bucket_sparse_sum += val1;
		}
		int newSample = discreteSample(bucket, K*2,
				rand01() * bucket_sparse_sum);
		if(newSample==-1){
			printf(" d=%d,l=%d | oldtopic=%d newtopic=%d\n",d,l ,k, newSample);
			for(int i=0;i<2*K;i++){
				printf("--bucket[%d]=%.14lf, fval=%.14lf\n",i,bucket[i],calculate_f_all_in_one(d,i/2,i%2==1));
			}
		}
		assert(newSample != -1);
		int newTopic=newSample/2;
		int newul=newSample%2;
		inc_mwk(w,newTopic);
		inc_ndk(d, newTopic);

		bool addRootTable=false,addTable=false;
//		addRootTable=(rand01()<=(1-(double)(t0k[newTopic])/(double)n0k[newTopic])?true:false);
		addRootTable=(addRootTable || t0k[newTopic]==0);
//		addTable=(rand01()<=(1-(double)(tdk[d][newTopic])/(double)ndk[d][newTopic])?true:false);
		addTable= (addTable || tdk[d][newTopic]==0);
		if(newul==0){
			if(addRootTable){inc_t0k(newTopic);}
			inc_n0k(newTopic);
			m0++;
		}
		if(addTable){
			inc_tdk(d, newTopic);
		}
		u[d][l] = newul;
		z[d][l] = newTopic;
	}
	virtual void gibbsStep() {
		for (size_t d = 0; d < docs.size(); d++) {
			VecInt& doc = docs[d];
			for (size_t l = 0; l < doc.size(); l++) {
				sampleWord(d, l);
			}
		}

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
		HDPLDA* lda=static_cast<HDPLDA*> (_lda);
		K=lda->K;
		V=lda->V;
		mwk=lda->mwk;
		mk=lda->mk;
		n0k=lda->n0k;
		t0k=lda->t0k;
		tk=lda->tk;
		tsum=lda->tsum;
		m0=lda->m0;
		b0=lda->b0;
		b1=lda->b1;
//		m0=0;
//		t0k.resize(K);
//		n0k.resize(K);
		gamma=lda->gamma;
		gammasum=lda->gammasum;
		gsnCache=lda->gsnCache;
		gsn_maxM=lda->gsn_maxM;
		gsn_maxN=lda->gsn_maxN;
	}
	virtual void copyState(void* _lda, std::string impl){
		if(impl=="HDPLDA"){
			HDPLDA* lda=static_cast<HDPLDA*> (_lda);
			K=lda->K;
			V=lda->V;
			mwk=lda->mwk;
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
		}else if(impl=="AliasHDPLDA"){
			AliasHDPLDA* lda=static_cast<AliasHDPLDA*> (_lda);
			K=lda->K;
			V=lda->V;
			mwk.resize(V);
			for(int w=0;w<V;w++){
				mwk[w].resize(K);
				for(int k=0;k<K;k++){
					mwk[w][k]=lda->getmwk(w,k);
				}
			}
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

	}
};
}

