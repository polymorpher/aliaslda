#pragma once
#include "GeneralizedStirlingNumber.hpp"
#include "ILDA.hpp"
#include"AliasWordPDPLDA.hpp"
namespace AliasLDA {
class WordPDPLDA: public ILDA {
	const int METRO_HASTING_STEPS = 2;
	const int ALIAS_SAMPLES_KEPT_RATIO = 8;
	const bool USE_MIXED_SUPPLY_CONSUME_MECH=true;

	Vector2DInt docs;
	int V;
	int K;
	double a; //discount
	double b; //concentration
	double gamma;  //for generating base Dirichlet
	double gammasum;
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
	Vector2DInt ndk;
	Vector2DInt z; //dim I _D __L. topic assignment
	VecInt nw;
	VecDouble bucket;
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;
	GeneralizedStirlingNumber& gsn;
	Vector2DDouble gsnCache;
	int gsn_maxN;
	int gsn_maxM;

	void initialiseGSN(){
		NMPair limits=gsn.getLimit(a);
		while(limits.m==-1)limits=gsn.getLimit(a);
		gsn_maxN=limits.n;
		gsn_maxM=limits.m;
		if(limits.m<gsn_maxM||limits.n<gsn_maxN){
			printf("initialising GSN Module...\n");
			while(!gsn.initialize(a,gsn_maxN,gsn_maxM)){}
			printf("GSN Module initialized\n");
		}
		double* tempCache;
		gsnCache.resize(gsn_maxN);
		while((tempCache=gsn.getCache(a))==NULL){}
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
		ndk[d][k]++;
	}
	inline void dec_ndk(int d, int k) {
		ndk[d][k]--;
	}

  	inline double getphi0(int w,int k){
  		return (gamma+t[w][k])/(gammasum+tsum[k]);
  	}
public:
	WordPDPLDA() :
			K(0), V(0),  alphasum(0), u01(0, 1), gsn(GeneralizedStirlingNumber::getInstance()) {
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
		ndk.resize(docs.size());
		for (int d = 0; d < docs.size(); d++) {
			ndk[d].resize(K);
		}
	}
	inline double calculate_f_flat(int w,int composite){
		return calculate_f(w,composite/2,composite%2==1);
	}
	inline double calculate_f(int w,int k, bool r){
		if(r){
//			printf("calc f true\n");
			return ((b+a*tsum[k])/(double)(b+msum[k]))*
						((t[w][k]+1)/(double)(m[w][k]+1))*
						getGSNRatio(m[w][k]+1,t[w][k]+1,m[w][k],t[w][k])*
						((gamma+t[w][k])/(double)(gammasum+tsum[k]));
		}else{
//			printf("calc f false\n");
			if(t[w][k]>0 && m[w][k]>0){
				return ((m[w][k]-t[w][k]+1)/(double)(m[w][k]+1))*getGSNRatio(m[w][k]+1,t[w][k],m[w][k],t[w][k])/(double)(b+msum[k]);

			}else{
				return 0;
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
	void initialiseGamma(){
		gammasum=gamma*V;
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
				if((t[w][k]==0) || (rand01() <= 1/(double)(t[w][k]+1))){
					inc_t(w,k);
				}
				assert(!(t[w][k]==0 && m[w][k]!=0));
			}
		}
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
		thetasum.resize(K);
		nw.resize(V);
		initialiseWordCount();
		initialiseDocumentTopicCounter();
		initialiseTopicsAndCounts();
		bucket.resize(2*K);
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
		double bucket_sparse_sum = 0;
		VecInt& ndkLocal = ndk[d];
		for (int i = 0; i < K; i++) {
			double val0 = (ndkLocal[i] +alpha) * calculate_f(w,i,false);
			bucket[2*i] = val0;
			bucket_sparse_sum += val0;

			double val1 = (ndkLocal[i] +alpha)  * calculate_f(w,i,true);
			bucket[2*i+1] = val1;
			bucket_sparse_sum += val1;
			assert(!(t[w][i]==0 && m[w][i]!=0));
		}
		int newSample = discreteSample(bucket, K*2,
				rand01() * bucket_sparse_sum);
		if(newSample==-1){
			printf(" d=%d,l=%d | oldtopic=%d newtopic=%d\n",d,l ,k, newSample);
			for(int i=0;i<2*K;i++){
				printf("--bucket[%d]=%lf, fval=%lf\n",i,bucket[i],calculate_f(w,i/2,i%2==1));
			}
		}
		assert(newSample != -1);
		int newTopic=newSample/2;
		int headOfTable=newSample%2;
		inc_m(w,newTopic);
		inc_ndk(d, newTopic);
		if(headOfTable==1){
			inc_t(w, newTopic);
		}
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
		WordPDPLDA* lda=static_cast<WordPDPLDA*> (_lda);
		K=lda->K;
		V=lda->V;
		m=lda->m;
		msum=lda->msum;
		alpha=lda->alpha;
		alphasum=lda->alphasum;
		alphaVec=lda->alphaVec;
		t=lda->t;
		tsum=lda->tsum;
		a=lda->a;
		b=lda->b;
		gamma=lda->gamma;
		gammasum=lda->gammasum;
		gsnCache=lda->gsnCache;
		gsn_maxM=lda->gsn_maxM;
		gsn_maxN=lda->gsn_maxN;
	}
	virtual void copyState(void* _lda, std::string impl){
		if(impl=="WordPDPLDA"){
			WordPDPLDA* lda=static_cast<WordPDPLDA*> (_lda);
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
		}else{
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
	}
};
}

